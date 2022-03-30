library(SingleCellExperiment)
library(doParallel)
library(N2R)
library(igraph)
library(RColorBrewer)


#List of functions to perform spatial sampling analysis 

#O)Auxiliary functions for color convertion :

color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("grey","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.99,na.rm = TRUE)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}


cluster_to_color = function(cluster_vector,Defined_list_cluster = NULL) {
  
  cluster_vector = as.character(cluster_vector)
  List_unique_cluster = unique(cluster_vector)
  List_unique_cluster = List_unique_cluster[order(List_unique_cluster)]
  
  if (!is.null(Defined_list_cluster)) {
    List_unique_cluster = Defined_list_cluster
  }
  
  N_clusters = length(unique(cluster_vector))
  
  optimal_palette = suppressWarnings(colorRamp(brewer.pal(N_clusters, "Spectral")))
  optimal_palette = optimal_palette((1:N_clusters)/N_clusters)
  optimal_palette = optimal_palette / 255
  optimal_palette = rgb(optimal_palette)
  
  color_cluster = cluster_vector
  
  for (k in 1:N_clusters) {
    selected_cluster = List_unique_cluster[k]
    color_cluster[cluster_vector==k] = optimal_palette[k]
  }
  return(color_cluster)
}


####Create_SCE() function :  Generates an SCE data container for multiplexed single-cell imaging data

Create_SCE = function (List_data, dimension = "2D", Bit_mode = 16, N_core = 6) 
{
  sce = SingleCellExperiment(assays = list(Raw_intensity = as.matrix(t(List_data$Expression_data))), 
                             colData = List_data$Cell_annotation, rowData = List_data$Gene_annotation, 
                             metadata = list(dimension = dimension, Bit_mode = Bit_mode, 
                                             N_core = N_core, Is_nuc_cyt = F))
  if ("Localisation" %in% colnames(rowData(sce))) {
    metadata(sce)$Is_nuc_cyt = T
  }
  return(sce)
}

###Count_normalization() function :  Normalisation and scaling of marker expression based on a glm count model
##sce	: a SingleCellExperiment.

#perform_batch_correction	: Boolean value, should covariates be regressed out ?
#batch_vector	: vector describing the batch of each cell. Required to perform the batch correction.
#residual_normalisation	: mathematical method for residual normalisation. Has to be chosen among "Anscombe","Pearson","Working" or "VST"

Count_normalization = function (sce, perform_batch_correction = FALSE, batch_vector = NULL, 
          residual_normalisation = "Anscombe") 
{
  if (!"Cell_size" %in% colnames(colData(sce))) {
    stop("The normalization procedure can not be performed as cell size is not available in the SCE object. Please select an other method or add a Cell_size column !")
  }
  if (!metadata(sce)$Is_nuc_cyt) {
    Cell_size = sce$Cell_size
    Transformed_data = t(assays(sce)[["Raw_intensity"]])
    Transformed_data = Transformed_data * 2^metadata(sce)$Bit_mode
    Transformed_data = apply(Transformed_data, MARGIN = 2, 
                             FUN = function(x) {
                               x * Cell_size
                             })
    Transformed_data = round(Transformed_data)
  }
  if (metadata(sce)$Is_nuc_cyt) {
    List_localisation = rowData(sce)[, "Localisation"]
    names(List_localisation) = rownames(sce)
    for (k in rownames(sce)) {
      if (List_localisation[k] == "Cytoplasm") {
        Object_size = sce$Cyto_size
      }
      if (List_localisation[k] == "Nuclear") {
        Object_size = sce$Nuc_size
      }
      if (List_localisation[k] == "Cell") {
        Object_size = sce$Cell_size
      }
      Transformed_data = t(assays(sce)[["Raw_intensity"]])
      Transformed_data = Transformed_data * 2^metadata(sce)$Bit_mode
      Transformed_data = apply(Transformed_data, MARGIN = 2, 
                               FUN = function(x) {
                                 x * Object_size
                               })
      Transformed_data = round(Transformed_data)
    }
  }
  cat(paste("Creating parallel backend using"), as.character(metadata(sce)$N_core), 
      "cores \n")
  registerDoParallel(metadata(sce)$N_core)
  if (!metadata(sce)$Is_nuc_cyt) {
    cat("Fitting Poisson regressions ...")
    List_regression_model = foreach(i = colnames(Transformed_data)) %dopar% 
      {
        Poisson_model = glm(Transformed_data[, i] ~ 
                              log(Cell_size), family = "poisson")
      }
  }
  cat(" done ! \n")
  if (metadata(sce)$Is_nuc_cyt) {
    cat("Fitting Poisson regressions ...")
    List_regression_model = foreach(k = colnames(Transformed_data)) %dopar% 
      {
        if (List_localisation[k] == "Cytoplasm") {
          Object_size = sce$Cyto_size
        }
        if (List_localisation[k] == "Nuclear") {
          Object_size = sce$Nuc_size
        }
        if (List_localisation[k] == "Cell") {
          Object_size = sce$Cell_size
        }
        Poisson_model = glm(Transformed_data[, k] ~ 
                              log(Object_size), family = "poisson")
      }
  }
  cat(" done ! \n")
  Residual_matrix = foreach(i = 1:length(List_regression_model), 
                            .combine = cbind) %dopar% {
                              Fitted_values = List_regression_model[[i]]$fitted.values
                              Real_values = Transformed_data[, i]
                              if (!residual_normalisation %in% c("Anscombe", "Pearson", 
                                                                 "Working", "VST", "Random_quantile")) {
                                cat("No proper method for residual normalization provided. Using the Anscombe normalization method")
                              }
                              if (residual_normalisation == "Anscombe") {
                                Normalised_residuals = 1.5 * (Real_values^(2/3) - 
                                                                Fitted_values^(2/3))/(Fitted_values^(1/6))
                              }
                              if (residual_normalisation == "Pearson") {
                                Normalised_residuals = (Real_values - Fitted_values)/sqrt(Fitted_values)
                              }
                              if (residual_normalisation == "Working") {
                                Normalised_residuals = (Real_values - Fitted_values)/Fitted_values
                              }
                              if (residual_normalisation == "VST") {
                                Normalised_residuals = (sqrt(Real_values) - sqrt(Fitted_values))/2
                              }
                              Normalised_residuals
                            }
  Residual_matrix = apply(Residual_matrix, MARGIN = 2, FUN = function(x) {
    x = x - min(x)
    x = x/max(x)
  })
  Residual_matrix = as.data.frame(Residual_matrix)
  colnames(Residual_matrix) = colnames(Transformed_data)
  assay(sce, "Count_normalised_intensity", withDimnames = FALSE) <- t(Residual_matrix)
  return(sce)
}


##### KNN_clustering() function : Perform cell clustering on a computed KNN graph

#sce	: a SingleCellExperiment.
#K : number of neighbors for the KNN graph computation
#clustering_method	: method used for graph clustering. Has to be chose among "Louvain","Greedy" and "Infomap"
#assay_type	: name of the data slot on which to perform the clustering (Raw_intensity, Arcsinh_transformed_intensity....). By default the regression-normalized data are used



KNN_clustering = function(sce, K = 30, clustering_method = "Louvain", assay_type = "Count_normalised_intensity", 
          metric = "angular") 
{
  if (!assay_type %in% names(assays(sce))) {
    stop("The slot required does not exist. Please select an existing slot !")
  }
  if (!clustering_method %in% c("Louvain", "Greedy", "Infomap")) {
    stop("The clustering method required does not exist. Please choose among Louvain,Greedy and Infomap !")
  }
  if (!metric %in% c("L2", "angular")) {
    stop("The distance metric required does not exist. Please choose among angular or L2!")
  }
  data_to_cluster = assay(sce, assay_type)
  data_to_cluster = t(data_to_cluster)
  Channel_for_clustering = rowData(sce)$Used_for_clustering
  if (sum(Channel_for_clustering) == 0) {
    Channel_for_clustering = rep(TRUE, ncol(data_to_cluster))
  }
  cat(paste("Clustering of the data using", as.character(sum(Channel_for_clustering)), 
            "channels from the", assay_type, "slot ! \n"))
  data_to_cluster = data_to_cluster[, Channel_for_clustering]
  KNN_graph_matrix = Knn(as.matrix(data_to_cluster), K, nThreads = metadata(sce)$N_core, 
                         verbose = F, indexType = metric)
  KNN_graph_matrix = KNN_graph_matrix + t(KNN_graph_matrix)
  Final_graph <- graph_from_adjacency_matrix(KNN_graph_matrix, 
                                             mode = "undirected", weighted = TRUE)
  cat("KNN computed \n")
  if (clustering_method == "Louvain") {
    Clustering <- cluster_louvain(Final_graph)
  }
  if (clustering_method == "Greedy") {
    Clustering <- cluster_fast_greedy(Final_graph, modularity = TRUE)
  }
  if (clustering_method == "Infomap") {
    Clustering <- cluster_infomap(Final_graph, modularity = TRUE)
  }
  Clustering_group = membership(Clustering)
  colLabels(sce) = Clustering_group
  cat(paste(as.character(length(unique(Clustering_group))), 
            "clusters have been identified \n"))
  return(sce)
}

### Plot_cluster_spatial() function : Function for visualisation the result of a clustering.
#sce	: a SingleCellExperiment. Clustering of the cell should have been previously performed.
#Image_number	: Number of name of the image/ROI to be plotted
#Cex_parameter :	Scaling factor for the size of the cells
#Specific_cluster	: Colors points based based on their belonging to a specific cluster
#Provided_group	: User-provided clustering

Plot_cluster_spatial = function (sce, Image_number = 1, Cex_parameter = 10, Specific_cluster = NULL, 
                                 Provided_cluster = NULL) 
{
  if (is.null(colLabels(sce))) {
    stop("Please compute clustering first ! \n")
  }
  Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber == 
                                                              Image_number], Y = sce$Location_Center_Y[sce$ImageNumber == 
                                                                                                         Image_number])
  Temp_cluster_data = colLabels(sce)[sce$ImageNumber == Image_number]
  if (!is.null(Provided_cluster)) {
    Temp_cluster_data = Provided_cluster
  }
  Dimension = metadata(sce)$dimension
  if ("Cell_size" %in% colnames(colData(sce))) {
    Temp_size_data = sce$Cell_size[sce$ImageNumber == Image_number]
  }
  if (!"Cell_size" %in% colnames(colData(sce))) {
    Temp_size_data = rep(10, sum(sce$ImageNumber == Image_number))
  }
  if (Dimension == "2D") {
    Temp_size_data = sqrt(Temp_size_data)
  }
  if (Dimension == "3D") {
    Temp_size_data = (Temp_size_data)^1/3
  }
  par(las = 1, bty = "l")
  if (is.null(Specific_cluster)) {
    color_temp_vector = cluster_to_color(Temp_cluster_data)
  }
  if (!is.null(Specific_cluster)) {
    if (Specific_cluster %in% unique(Temp_cluster_data)) {
      color_temp_vector = (as.numeric(Temp_cluster_data == 
                                        Specific_cluster))
    }
  }
  plot(Temp_location_data, cex = Temp_size_data/Cex_parameter, 
       pch = 21, bg = color_temp_vector)
}

##Plot_gene_expression_spatial() function :Plot the spatial distribution of a given gene normalised expression

#sce: a SingleCellExperiment object.
#Image_number	: Number of name of the image/ROI to be plotted
#Cex_parameter	: Scaling factor for the size of the cells
#assay_type	:name of the data slot on which to perform the clustering (Raw_intensity, Arcsinh_transformed_intensity....). By default the regression-normalized data are used
#Gene	: Name of the gene to be plotted


Plot_gene_expression_spatial = function (sce, Image_number = 1, Cex_parameter = 5, assay_type = "Count_normalised_intensity", 
          Gene = NULL) 
{
  if (is.null(Gene) | !Gene %in% rownames(sce@assays@data@listData$Raw_intensity)) {
    stop("Please select a correct gene to plot ! \n")
  }
  if (!assay_type %in% names(assays(sce))) {
    stop("The slot required does not exist. Please select an existing slot !")
  }
  Dimension = metadata(sce)$dimension
  if ("Cell_size" %in% colnames(colData(sce))) {
    Temp_size_data = sce$Cell_size[sce$ImageNumber == Image_number]
  }
  if (!"Cell_size" %in% colnames(colData(sce))) {
    Temp_size_data = rep(10, sum(sce$ImageNumber == Image_number))
  }
  if (Dimension == "2D") {
    Temp_size_data = sqrt(Temp_size_data)
  }
  if (Dimension == "3D") {
    Temp_size_data = (Temp_size_data)^1/3
  }
  Temp_location_data = data.frame(X = sce$Location_Center_X[sce$ImageNumber == 
                                                              Image_number], Y = sce$Location_Center_Y[sce$ImageNumber == 
                                                                                                         Image_number])
  Temp_expression_data = as.numeric(assay(sce, assay_type)[Gene, 
  ])
  Temp_expression_data = Temp_expression_data - min(Temp_expression_data)
  par(las = 1, bty = "l")
  plot(Temp_location_data, cex = Temp_size_data/Cex_parameter, 
       pch = 21, bg = color_convertion(Temp_expression_data)[sce$ImageNumber == 
                                                                Image_number])
}


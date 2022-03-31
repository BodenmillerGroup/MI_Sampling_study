#List of functions to perform spatial sampling analysis 

library("SingleCellExperiment")
library("doParallel")
library("RColorBrewer")
library("CountClust")
library("N2R")
library("igraph")
library("SQUAREM")


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



##I)

#II)Functions for the sanmpling analysis in itself 

#Random_spatial_sampling() function : basic function that randomly draw N_samplings rectangular ROIs of 
#width width_FOV and of height height_FOV from a large spatial experiment 


#sce:	a SingleCellExperiment.
#width_FOV: height of the individual rectangles
#height_FOV	: width of the individual rectangles
#N_samplings: number of individual rectangles samples
#Selected_image	: of which image/FOV should the sampling be performed ?
  

Random_spatial_sampling = function(sce, width_FOV = 400, height_FOV = 400, N_samplings = 10, 
          Selected_image = 1, plot_result = TRUE) 
{
  sce = sce[, sce$ImageNumber == Selected_image]
  Equivalent_radius = sqrt((width_FOV^2)/4 + (height_FOV^2)/4) * 
    2
  x_range_sampling = range(sce$Location_Center_X)
  y_range_sampling = range(sce$Location_Center_Y)
  List_center = c()
  for (k in 1:N_samplings) {
    Is_in_empty_space = FALSE
    while (!Is_in_empty_space) {
      position_temp = c(runif(n = 1, min = x_range_sampling[1], 
                              max = x_range_sampling[2]), runif(n = 1, min = y_range_sampling[1], 
                                                                max = y_range_sampling[2]))
      Dist_matrix = dist(rbind(position_temp, List_center), 
                         method = "manhattan")
      Dist_matrix = as.matrix(Dist_matrix)
      if (nrow(Dist_matrix) == 1) {
        Is_in_empty_space = T
        List_center = rbind(List_center, position_temp)
      }
      if (nrow(Dist_matrix) != 1) {
        List_distance_temp = Dist_matrix[1, ]
        if (min(List_distance_temp[-1]) > Equivalent_radius) {
          Is_in_empty_space = T
          List_center = rbind(List_center, position_temp)
        }
      }
    }
  }
  List_sampled_cells = c()
  List_sampled_cluster = c()
  List_sample = c()
  for (k in 1:nrow(List_center)) {
    center_temp = List_center[k, ]
    Selected_cells = which(sce$Location_Center_X > (center_temp[1] - 
                                                      width_FOV/2) & sce$Location_Center_X < center_temp[1] + 
                             width_FOV/2 & sce$Location_Center_Y > center_temp[2] - 
                             height_FOV/2 & sce$Location_Center_Y < center_temp[2] + 
                             height_FOV/2)
    List_sampled_cells[[k]] = Selected_cells
    List_sampled_cluster[[k]] = colLabels(sce)[Selected_cells]
    List_sample = c(List_sample, rep(paste("Sample", k, 
                                           sep = "_"), length(Selected_cells)))
  }
  if (plot_result) {
    par(bty = "n", las = 1)
    plot(sce$Location_Center_X, sce$Location_Center_Y, pch = 21, 
         bg = cluster_to_color(colLabels(sce)))
    for (k in 1:nrow(List_center)) {
      center_temp = List_center[k, ]
      rect(xleft = center_temp[1] - width_FOV/2, ybottom = center_temp[2] - 
             height_FOV/2, xright = center_temp[1] + width_FOV/2, 
           ytop = center_temp[2] + height_FOV/2, col = "black", 
           density = 40)
    }
  }
  List_result = list(List_sampling = List_sample, List_sampled_cells = unlist(List_sampled_cells), 
                     List_sampled_cluster = unlist(List_sampled_cluster))
  return(List_result)
}

#####

###Perform_sampling_analysis() function : Perform a sampling analysis where one or various sampling parameters are varying.
#

#sce :a SingleCellExperiment.
#Selected_image	:of which image/FOV should the sampling analysis be performed ?
#N_times : number of times each type of sampling is performed (typically between 20 and 100)
#N_sampling_region_vector	 : vector (or integer) describing the different values taken by the number of sampled regions
#width_FOV_vector	 : vector (or real number) describing the width of the FOV
#height_FOV_vector : vector (or real number) describing the height of the FOV
#Threshold_detection_cluster: real number corresponding to the minimal number of a given cell type to be considered as detected. Typically around 50 or 100 for IMC data, rather 2/3 for Visium


Perform_sampling_analysis = function (sce, Selected_image = 1, N_times = 50, N_sampling_region_vector = 1:20, 
          width_FOV_vector = 400, height_FOV_vector = 400, Threshold_detection_cluster = 50) 
{
  Length_vector = c(length(N_sampling_region_vector), length(width_FOV_vector), 
                    length(height_FOV_vector))
  Length_vector = unique(Length_vector)
  Length_vector = Length_vector[Length_vector != 1]
  if (length(Length_vector) > 1) {
    stop("multiple parameter vectors of with a size bigger than 1 have been provided. Please use appropriate parameters !")
  }
  if (length(N_sampling_region_vector) == 1) {
    N_sampling_region_vector = rep(N_sampling_region_vector, 
                                   Length_vector)
  }
  if (length(width_FOV_vector) == 1) {
    width_FOV_vector = rep(width_FOV_vector, Length_vector)
  }
  if (length(height_FOV_vector) == 1) {
    height_FOV_vector = rep(height_FOV_vector, Length_vector)
  }
  Global_composition = table(factor(colLabels(sce)))
  Global_composition_normalised = Global_composition/sum(Global_composition)
  Global_composition = log(Global_composition/prod(Global_composition)^(1/length(Global_composition)))
  Mean_number_cluster_identified = c()
  Sd_number_cluster_identified = c()
  Mean_divergence_global_composition = c()
  Sd_divergence_global_composition = c()
  Mean_correlation_composition = c()
  Sd_correlation_composition = c()
  Mean_KL_score = c()
  Sd_KL_score = c()
  for (i in 1:Length_vector) {
    print(i)
    Number_cluster_identified_temp = c()
    Divergence_temp = c()
    Correlation_temp = c()
    KL_temp = c()
    for (j in 1:N_times) {
      x = Random_spatial_sampling(sce, Selected_image = Selected_image, 
                                  width_FOV = width_FOV_vector[i], height_FOV = height_FOV_vector[i], 
                                  N_samplings = N_sampling_region_vector[i], plot_result = F)
      table_sampled_clusters = table(factor(x$List_sampled_cluster, 
                                            levels = levels(factor(colLabels(sce)))))
      table_sampled_clusters_normalized = table_sampled_clusters/sum(table_sampled_clusters)
      Number_cluster_identified_temp = c(Number_cluster_identified_temp, 
                                         sum(table_sampled_clusters > Threshold_detection_cluster))
      KL_temp = c(KL_temp, sum(table_sampled_clusters_normalized * 
                                 log(table_sampled_clusters_normalized/Global_composition_normalised), 
                               na.rm = T))
      sampled_composition = table_sampled_clusters + 1
      sampled_composition = log(sampled_composition/(prod(sampled_composition)^(1/length(sampled_composition))))
      Aitchison_distance = sqrt(sum((sampled_composition - 
                                       Global_composition)^2))
      Divergence_temp = c(Divergence_temp, Aitchison_distance)
      table_sampled_cluster_normalised = table_sampled_clusters/sum(table_sampled_clusters)
      R_temp = cor(table_sampled_cluster_normalised, Global_composition_normalised)
      Correlation_temp = c(Correlation_temp, R_temp)
    }
    Mean_divergence_global_composition = c(Mean_divergence_global_composition, 
                                           mean(Divergence_temp))
    Sd_divergence_global_composition = c(Sd_divergence_global_composition, 
                                         sd(Divergence_temp))
    Mean_number_cluster_identified = c(Mean_number_cluster_identified, 
                                       mean(Number_cluster_identified_temp))
    Sd_number_cluster_identified = c(Sd_number_cluster_identified, 
                                     sd(Number_cluster_identified_temp))
    Mean_correlation_composition = c(Mean_correlation_composition, 
                                     mean(Correlation_temp))
    Sd_correlation_composition = c(Sd_correlation_composition, 
                                   sd(Correlation_temp))
    Mean_KL_score = c(Mean_KL_score, mean(KL_temp))
    Sd_KL_score = c(Sd_KL_score, sd(KL_temp))
  }
  Statistic_data_frame = data.frame(N_sampling = N_sampling_region_vector, 
                                    Mean_number_cluster = Mean_number_cluster_identified, 
                                    Sd_number_cluster = Sd_number_cluster_identified, Mean_KL_divergence = Mean_KL_score, 
                                    Sd_KL_divergence = Sd_KL_score, Mean_correlation_composition = Mean_correlation_composition, 
                                    Sd_correlation_composition = Sd_correlation_composition, 
                                    Mean_Aitchison_distance = Mean_divergence_global_composition, 
                                    Sd_Aitchison_distance = Sd_divergence_global_composition)
  return(Statistic_data_frame)
}

#####

##Visualize_simple_sampling() function: Function that plot and analyze the results from the Perform_sampling_analysis() and estimate the tau and No parameters.

#Sampling_data_frame : a data.frame object resulting from the Perform_sampling_analysis() function

Visualize_simple_sampling = function (Sampling_data_frame) 
{
  Max_cluster_recovered = round(1.2 * max(Sampling_data_frame$Mean_number_cluster))
  Max_sampled_FoV = round(1.2 * max(Sampling_data_frame$N_sampling))
  par(las = 1, bty = "l")
  plot(Sampling_data_frame$N_sampling, Sampling_data_frame$Mean_number_cluster, 
       xlim = c(0, Max_sampled_FoV), ylim = c(0, Max_cluster_recovered), 
       xaxs = "i", yaxs = "i", xlab = "Number of sampled regions", 
       ylab = "Mean number of recovered clusters", cex.lab = 1.2, 
       pch = 21, bg = "red3", cex = 2)
  y = Sampling_data_frame$Mean_number_cluster
  x = Sampling_data_frame$N_sampling
  expo_model_number_cluster = nls(y ~ N * (1 - exp(-x/tau)), 
                                  start = list(N = 20, tau = 5))
  curve(expr = coef(expo_model_number_cluster)[1] * (1 - exp(-x/coef(expo_model_number_cluster)[2])), 
        add = T, col = "red", lwd = 2, from = 0, to = 100, lty = 2)
  Result_vector = coef(expo_model_number_cluster)
  R_squared = cor(predict(expo_model_number_cluster, newdata = x), 
                  y)^2
  Result_vector = c(Result_vector, R_squared)
  names(Result_vector)[3] = "R_squared"
  abline(h = Result_vector[1], lwd = 1.5, lty = 2, col = "grey")
  legend("right", legend = c(paste("N =", round(Result_vector[1], 
                                                2)), paste("tau =", round(Result_vector[2], 2)), paste("R2 =", 
                                                                                                       round(Result_vector[3], 3))), bty = "n", cex = 1.2)
  return(Result_vector)
}


#####

#Visualize_complex_sampling() function : Function that plot and analyze the results from the Perform_sampling_analysis() and estimate the tau and No parameters when various parameters value have been used .

#Sampling_data_frame	: a data.frame object resulting from the Perform_sampling_analysis() function
#Parameter_table : a table that describe the different parameter values


Visualize_complex_sampling = function (Sampling_data_frame, Parameter_table) 
{
  Max_cluster_recovered = round(1.2 * max(Sampling_data_frame$Mean_number_cluster))
  Max_sampled_FoV = round(1.2 * max(Sampling_data_frame$N_sampling))
  X = c()
  for (k in 1:ncol(Parameter_table)) {
    X = paste(X, Parameter_table[, k], sep = "_")
  }
  X = substr(X, start = 2, stop = 500)
  Parameter_values = unique(X)
  par(las = 1, bty = "l")
  plot(Sampling_data_frame$N_sampling, Sampling_data_frame$Mean_number_cluster, 
       xlim = c(0, Max_sampled_FoV), ylim = c(0, Max_cluster_recovered), 
       xaxs = "i", yaxs = "i", xlab = "Number of sampled regions", 
       ylab = "Mean number of recovered clusters", cex.lab = 1.2, 
       pch = 21, bg = cluster_to_color(as.numeric(as.factor(X))), 
       cex = 2)
  segments(x0 = Sampling_data_frame$N_sampling, x1 = Sampling_data_frame$N_sampling, 
           y0 = Sampling_data_frame$Mean_number_cluster - Sampling_data_frame$Sd_number_cluster, 
           y1 = Sampling_data_frame$Mean_number_cluster + Sampling_data_frame$Sd_number_cluster, 
           lwd = 0.2)
  Result_table = c()
  for (k in 1:length(Parameter_values)) {
    y = Sampling_data_frame$Mean_number_cluster[X == Parameter_values[k]]
    x = Sampling_data_frame$N_sampling[X == Parameter_values[k]]
    expo_model_number_cluster = nls(y ~ N * (1 - exp(-x/tau)), 
                                    start = list(N = 20, tau = 5))
    curve(expr = coef(expo_model_number_cluster)[1] * (1 - 
                                                         exp(-x/coef(expo_model_number_cluster)[2])), add = T, 
          col = "red", lwd = 2, from = 0, to = 100, lty = 2)
    Result_vector = coef(expo_model_number_cluster)
    R_squared = cor(predict(expo_model_number_cluster, newdata = x), 
                    y)^2
    Result_vector = c(Result_vector, R_squared)
    names(Result_vector)[3] = "R_squared"
    Result_table = rbind(Result_table, Result_vector)
  }
  abline(h = max(Result_table[, 1]), lwd = 2, lty = 2, col = "grey")
  return(Result_table)
}


###

# Global_alpha_estimation : Function that estimates the alpha parameter of a given tissue using many small FoVs instead of a unique large FoVs.

Global_alpha_estimation = function (sce, List_image = NULL, Threshold_detection = 20) 
{
  if (is.null(Value_split)) {
    Value_split = c(1, 1.1, 1.2, 1.5, 1.8, 2, 2.5, 3)
  }
  if (is.null(List_image)) {
    List_image = unique(sce$ImageNb)
  }
  List_coef = c()
  R_squared = c()
  for (k in List_image) {
    sce_temp = sce[, sce$ImageNb == k]
    x_range = max(sce_temp$Location_Center_X)
    y_range = max(sce_temp$Location_Center_Y)
    center_x = x_range/2
    center_y = y_range/2
    N_clusters_detected = c()
    for (i in Value_split) {
      x_range_temp = c(center_x - x_range/(i * 2), center_x + 
                         x_range/(i * 2))
      y_range_temp = c(center_y - y_range/(i * 2), center_y + 
                         y_range/(i * 2))
      sce_temp_filtered = sce_temp[, sce_temp$Location_Center_X < 
                                     x_range_temp[2] & sce_temp$Location_Center_X > 
                                     x_range_temp[1] & sce_temp$Location_Center_Y > 
                                     y_range_temp[1] & sce_temp$Location_Center_Y < 
                                     y_range_temp[2]]
      x = table(factor(colLabels(sce_temp_filtered), levels = unique(colLabels(sce))))
      N_clusters_detected = c(N_clusters_detected, sum(x > 
                                                         Threshold_detection))
    }
    m = lm(log10(N_clusters_detected) ~ log10(Value_split), 
           subset = N_clusters_detected > 0)
    List_coef = c(List_coef, coef(m)[2])
    R_squared = c(R_squared, summary(m)$r.squared)
  }
  Table_estimation = data.frame(Alpha = -List_coef, R_squared = R_squared, 
                                row.names = List_image)
  return(Table_estimation)
}









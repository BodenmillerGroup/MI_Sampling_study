library(pagoda2)
library(CountClust)
library(doParallel)
library(SingleCellExperiment)

#O)Auxiliary functions for visualisation :

string.to.colors = function (string, colors = NULL) 
{
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      (break)("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  }
  else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x) {
    conv[which(conv[, 1] == x), 2]
  }))
}



color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("white","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.999,na.rm = T)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}




#I) Data downloading, loading and processing 

#A)Downloading the raw file 

#Here as an example :

setwd("~/Desktop/Scripts_publication_NM/Test/")

#Downloading the expression matrix : 
download.file(url = "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_raw_feature_bc_matrix.tar.gz",destfile = "V1_Human_Lymph_Node_raw_feature_bc_matrix.tar.gz")
system("tar -xvf V1_Human_Lymph_Node_raw_feature_bc_matrix.tar.gz")
system("rm V1_Human_Lymph_Node_raw_feature_bc_matrix.tar.gz")

#Downloading the spatial metadata files (H&E images and barcodes/location correspondance): 
download.file(url = "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_spatial.tar.gz",destfile = "V1_Human_Lymph_Node_spatial.tar")
system("tar -xvf V1_Human_Lymph_Node_spatial.tar")
system("rm V1_Human_Lymph_Node_spatial.tar")


#B)Loading the data and metadata

data_raw = read.10x.matrices("raw_feature_bc_matrix/")
l = colnames(data_raw)
l = substr(l,start = 5,stop = 600)
colnames(data_raw) = l

data_location = read.delim("spatial//tissue_positions_list.csv",sep=",",header = F,row.names = 1)
data_location = data_location[colnames(data_raw),]
data_location = data_location[,c(2,3)]
colnames(data_location) = c("Location_X","Location_Y")
data_location$Location_X = data_location$Location_X*6400/78
data_location$Location_Y = data_location$Location_Y*6400/128

#C)Filtering

Lib_size = colSums(data_raw)
hist(log10(Lib_size),100)
Gene_size = rowSums(data_raw)
hist(log10(Gene_size+1),100)

par(las=1,bty="l")
plot(data_location,pch=21,bg=string.to.colors((Lib_size<1000)))

data_count = data_raw[Gene_size>10,Lib_size>1000]

Multiple_gene = names(which(table(rownames(data_count))>1))
data_count = data_count[!rownames(data_count)%in%Multiple_gene,]
Location_count = data_location[Lib_size>1000,]

#II) Data analysis

#A)Creating the Pagoda2 object and the LDA dimensionality reduction
r <- Pagoda2$new(data_count,log.scale=T)
r$adjustVariance(plot=T,gam.k=5)

#B)Latent Dirichlet Analysis (LDA)

#Selecting the top 1500 genes and removing 'bad genes' : Immunoglobulin and Mitochondrial genes
Selected_genes = r$getOdGenes(1500)
Selected_genes = Selected_genes[!grepl(Selected_genes,pattern = "IGH")]
Selected_genes = Selected_genes[!grepl(Selected_genes,pattern = "IGK")]
Selected_genes = Selected_genes[!grepl(Selected_genes,pattern = "IGL")]
Selected_genes = Selected_genes[!grepl(Selected_genes,pattern = "MT-")]


#Computational intensive part : computing LDA with 5 different number of topics in a parallel fashion to know which number of topic should be selected 
registerDoParallel(cores=5)

Model_LDA_merged = foreach(k=c(5,10,15,20,25)) %dopar% {
  Model_LDA_temp = FitGoM(as.matrix(t(data_count[Selected_genes,])),K = k,tol=100,options="BIC")
  Model_LDA_temp
}


List_BIC = unlist(lapply(Model_LDA_merged,FUN = function(x) {x$BIC}))

#How many topics should we select ? 
plot(c(5,10,15,20,25),List_BIC/10000,pch=21,bg="orange",type="o",cex=2,xlab="Number of topics",ylab="BIC",cex.lab=1.3,cex.axis=1)

#Here : 15 topics is enough 
Selected_models = Model_LDA_merged[[3]]
Mixing = Selected_models$fit$omega
Contribution = Selected_models$fit$theta
Marginal_distribution = rowSums(data_count[Selected_genes,])/sum(data_count[Selected_genes,])

Get_relevance_table = function(Contribution,Marginal_distribution,lambda=0) {
  Relevance_table = apply(Contribution,MARGIN = 2,FUN = function(x) {lambda*log(x) + (1-lambda)*log(x/Marginal_distribution)})
}

#This table allows to biologically interpret each 
Contribution_normalized = Get_relevance_table(Contribution,Marginal_distribution,lambda=0)


#C) Performing clustering on the LDA result

r$reductions$LDA = Mixing
r$makeKnnGraph(k=15,type='LDA',distance = "cosine")
r$getKnnClusters(method=multilevel.community,type='LDA',name = "LDA_cluster")
r$getDifferentialGenes(type = "LDA",clusterType = "LDA_cluster",verbose = T,z.threshold = 3)

plot(Location_count,pch=16,col=string.to.colors(r$clusters$LDA$LDA_cluster))

#III) Creating an SCE object

sce = SingleCellExperiment(assays = list(Raw_intensity = as.matrix(t(Mixing))))
colLabels(sce) = as.numeric(r$clusters$LDA$LDA_cluster)
sce$Location_Center_X = Location_count$Location_X
sce$Location_Center_Y = Location_count$Location_Y
sce$ImageNumber =  1

#IV)Testing the spatial sampling functions


Random_spatial_sampling(sce,width_FOV = 400,height_FOV = 400,N_samplings = 10)
Simple_sampling_analysis = Perform_sampling_analysis(sce,Selected_image = 1,N_times = 50,N_sampling_region_vector = 1:20,
                                                     width_FOV_vector = 400,height_FOV_vector = 400,Threshold_detection = 2)

Visualize_simple_sampling(Simple_sampling_analysis)


height_vector = rep(c(200,300,400,500,600),each=30)
width_vector = rep(c(200,300,400,500,600),each=30)
N_sampling_region_vector = rep((1:30),5)

Complex_sampling = Perform_sampling_analysis(sce,Selected_image = 1,N_times = 50,
                                             N_sampling_region_vector = N_sampling_region_vector,
                                             width_FOV_vector = width_vector,
                                             height_FOV_vector = height_vector,
                                             Threshold_detection = 2)


Parameter_table = data.frame(Height =height_vector,
                             Width =width_vector)

Fitting_tau = Visualize_complex_sampling(Complex_sampling,Parameter_table)


FoV_width = c(200,300,400,500,600)
plot(FoV_width,Fitting_tau[,"tau"],log="xy",xlim=c(150,700),ylim=c(1,30),
     xaxs='i',yaxs='i',xlab="FoV's width",ylab="Tau parameter",cex.lab=1.3,pch=21,bg="red3",cex=2)
m = lm(log10(Fitting_tau[,"tau"])~log10(FoV_width))
abline(coef(m),lwd=2,lty=2,col="grey")
print(coef(m))



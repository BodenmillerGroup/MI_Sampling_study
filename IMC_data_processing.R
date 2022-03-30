

source("Desktop/Scripts_publication_NM/List_scripts_IMC_processing.R")

Cell_information = read.delim("Steinbock_directory/regionprops/20210805_LN_panorama_3_1.csv",sep=",")
Cell_information_reshaped = data.frame(ImageNumber = 1,ObjectNumber = Cell_information$Object,
                                       Location_Center_X = Cell_information$centroid.0,Location_Center_Y = Cell_information$centroid.1,
                                       Cell_size = Cell_information$area)
Protein_expression = read.delim("Steinbock_directory/intensities/20210805_LN_panorama_3_1.csv",sep=",")

Protein_expression = Protein_expression[,-1]

#B)Merging
Gene_annotation = data.frame(Used_for_clustering = rep(FALSE,ncol(Protein_expression)),
                             DNA_channel =  rep(FALSE,ncol(Protein_expression)),row.names =colnames(Protein_expression) )

Gene_annotation$Used_for_clustering = !rownames(Gene_annotation)%in%c("DNA_1","DNA_2","Ki.67","Cleaved.Caspase3")



List_data_panorama = list(Expression_data = Protein_expression,
                          Cell_annotation = Cell_information_reshaped,
                          Gene_annotation = Gene_annotation)
sce = Create_SCE(List_data_panorama,dimension = "2D",Bit_mode = 32,N_core = 8)


sce = Count_normalization(sce,residual_normalisation = "Anscombe")
sce = KNN_clustering(sce,K = 15,clustering_method = "Louvain",assay_type = "Count_normalised_intensity",metric = "L2")

Plot_cluster_spatial(sce=sce,Image_number = 1,Cex_parameter = 12)
Plot_gene_expression_spatial(sce=sce,Image_number = 1,Cex_parameter = 12,Gene = "FOXP3")

###

source("Desktop/Scripts_publication_NM/List_scripts_sampling.R")
Simple_sampling_analysis = Perform_sampling_analysis(sce,Selected_image = 1,N_times = 50,N_sampling_region_vector = 1:20,width_FOV_vector = 400,height_FOV_vector = 400,Threshold_detection = 50)
Visualize_simple_sampling(Simple_sampling_analysis)



height_vector = rep(c(200,250,300,350,400,450,500),each=10)
width_vector = rep(c(200,250,300,350,400,450,500),each=10)
N_sampling_region_vector = rep((1:10),7)

Complex_sampling = Perform_sampling_analysis(sce,Selected_image = 1,N_times = 50,
                                             N_sampling_region_vector = N_sampling_region_vector,
                                             width_FOV_vector = width_vector,
                                             height_FOV_vector = height_vector,
                                             Threshold_detection = 50)


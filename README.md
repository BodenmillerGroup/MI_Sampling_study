# Spatial sampling analysis of Multiplexed Imaging (MI) and Spatial Transcriptomic (ST) data

This repository contains all the code and scripts used to perform all the spatial sampling analysis described in ["Optimizing multiplexed imaging experimental design through tissue spatial segregation estimation"](https://www.biorxiv.org/content/10.1101/2021.11.28.470262v2). On a practical level, the scripts rely on the [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) objects and can be used on both ST (i.e Visium) and MI (IMC, MIBI or CODEX) datasets and only require minimal pre-processing.

This repository contains the following files :

- **List_scripts_sampling.R** : file containing all the functions required for the spatial sampling analysis.
- **Visium_data_processing.R** : script download, process and analyse Visium® datasets.
- **Example_sce_IMC.rds** : a rds file that contains a processed sce object derived from the IMC analysis of a healthy human lymph node. In order to have a lightweight file, only the cell type information is available for each cell, not the individual protein expression.
-  **Example_sce_Visium.rds** : a rds file that contains a processed sce object derived from the Visium analysis of a human healthy human lymph node. In order to have a lightweight file, only the cell type information is available for each cell, the individual spot RNA expression profile is not provided.


## Installing required packages 

Several packages are needed to perform the different scripts. They can be installed by the following command :

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment","doParallel","RColorBrewer","CountClust","N2R","igraph"))
```

In addition, the [Pagoda2 package](https://github.com/kharchenkolab/pagoda2) has to be installed to process the spatial transcriptomic data in an efficient manner. Please check on the corresponding github page the different dependencies needed to install it. 


## Pre-processing the data 

Data has to be stored in a SingleCellExperiment (SCE) object in order to be properly analysed. For a thoroughly description and introduction to the SCE object please have a look at this [introduction](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html). The following constraints need to be fullfill :

- The cell/spot classification/clustering need to be stored in the **ColLabels** part of the SCE object. It has to be a numerical vector for compatibility.
- X and Y location of cells should be stored in the column metadata (ColData) of the SCE object with the **"Location_Center_X"** and **"Location_Center_Y"** names. In addition the **"ImageNumber"** column will be added and set to 1 (all cells are supposed to come from the same image).

Scripts to transform raw Visium or MI datasets into a usable SCE object are available in this repository (see above).

## Performing regular sampling analysis using one large Field of View (FoV)

We assume a properly organised SCE object (called here sce) is available. In addition the **List_scripts_sampling.R** file has been downloaded locally.
We start by loading the different functions from the R file :

```r
source("path/to/file/List_scripts_sampling.R")
```


We can then perform a basic sampling analysis :

```r
Simple_sampling_analysis = Perform_sampling_analysis(sce,Selected_image = 1,N_times = 50,N_sampling_region_vector = 1:20,width_FOV_vector = 400,height_FOV_vector = 400,Threshold_detection = 50)
```
Here various samplings are tested with different number of ROIs (from 1 to 20) and with ROIs being squares of 400µm. We also specify that 50 cells/spots of each group/cluster have to be sampled in order to consider a group to be detected.

We can then fit the empirical model described in the manuscript and extract the two model parameters :

```r
Visualize_simple_sampling(Simple_sampling_analysis)
```

More complex analyses can also be performed. For instance we can measure the impact of the ROIs size on the tau parameter. We thus perform sampling with 7 different values 

```r
height_vector = rep(c(200,250,300,350,400,450,500),each=10)
width_vector = rep(c(200,250,300,350,400,450,500),each=10)
N_sampling_region_vector = rep((1:10),7)

Complex_sampling = Perform_sampling_analysis(sce,Selected_image = 1,N_times = 50,
                                             N_sampling_region_vector = N_sampling_region_vector,
                                             width_FOV_vector = width_vector,
                                             height_FOV_vector = height_vector,
                                             Threshold_detection = 50)

```
This step can be quite computationally heavy and will likely take some time. Once it is finished you can extract the different tau values easily :

```r
Parameter_table = data.frame(Height =height_vector,
                             Width =width_vector)

Fitting_tau = Visualize_complex_sampling(Complex_sampling,Parameter_table)
```

## Performing sampling analysis using a large number of small FoVs

To conclude this tutorial, we will see how to estimate the 𝛂 parameter using a large number of small FoV. This can be done using large-scale MI datasets for instance. We assume here that a correct SCE object is available and that it contains data from multiple FoVs. The FoV will be encoded in the **ImageNumber** channel (numerical values).

First we obtain a rough alpha estimate for each FoV :

```r
Alpha_estimate = Global_alpha_estimation(sce)
```
This generates a data.table object where the first column corresponds to the alpha estimate and the second to the quality of the estimate (correlation coefficient). Estimations with a low coefficients correlation (i.e unreliable/low-quality) should be removed and only high quality estimates used for final estimation :

```r
mean(Alpha_estimate$Alpha[Alpha_estimate$R_squared>0.9])
```


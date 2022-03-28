# Spatial sampling analysis of Multiplexed Imaging (MI) and Spatial Transcriptomic (ST) data

This repository contains all the code and scripts used to perform all the spatial sampling analysis described in ["Optimizing multiplexed imaging experimental design through tissue spatial segregation estimation"](https://www.biorxiv.org/content/10.1101/2021.11.28.470262v2). On a practical level, the scripts rely on the [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) objects and can be used on both ST (i.e Visium) and MI (IMC, MIBI or CODEX) datasets and only require minimal pre-processing.

This repository contains the following scripts :

- **List_scripts_sampling.R** : file containing all the functions required for the spatial sampling analysis
- **Visium_data_processing.R** : example of how to download, process and analyse Visium® datasets.

## Installing required packages 

Several packages are needed to perform the different scripts. They can be installed by the following command :

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment","doParallel","RColorBrewer","CountClust))
```

In addition, the [Pagoda2 pacakge](https://github.com/kharchenkolab/pagoda2) has to be installed to proces. Please check on the corresponding github page the different dependencies needed to install it. 


## Pre-processing the data 

Data has to be stored in a SingleCellExperiment (SCE) object in order to be properly analysed. For a thoroughly description and introduction to the SCE object please have a look at this [introduction](https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html). The following constraints need to be fullfill :

- The cell/spot classification/clustering need to be stored in the ColLabels part of the SCE object. It has to be a numerical vector for compatibility.
- X and Y location of cells should be stored in the column metadata (ColData) of the SCE object with the "Location_Center_X" and "Location_Center_Y" names. In addition the "ImageNumber" column will be added and set to 1 (all cells are supposed to come from the same image).

Scripts to transform raw Visium or MI datasets into a usable SCE object are available in this repository (see above).

## Performing sampling analysis

We assume a properly organised SCE object (called here sce) is available. In addition the **List_scripts_sampling.R** file has been downloaded locally.
We start by loading the different functions from the R file :

```r
source("path/to/file/List_scripts_sampling.R")
```

We can then perform 






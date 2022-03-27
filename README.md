# Spatial sampling analysis of Multiplexed Imaging (MI) and Spatial Transcriptomic (ST) data

This repository contains all the code and scripts used to perform all the spatial sampling analysis described in ["Optimizing multiplexed imaging experimental design through tissue spatial segregation estimation"](https://www.biorxiv.org/content/10.1101/2021.11.28.470262v2). On a practical level, the scripts rely on the [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html) objects and can be used on both ST (i.e Visium) and MI (IMC, MIBI or CODEX) datasets and only require minimal pre-processing.

This repository contains the following scripts :

- **List_scripts_sampling.R** : file containing all the functions required for the spatial sampling analysis
- **Visium_data_processing.R** : example of how to download, process and analyse Visium® datasets.

## Installing required packages 

Several packages are needed to perform the different scripts. They can be installed by the following command :

```r
BiocManager::install(c("SingleCellExperiment","doParallel","RColorBrewer","CountClust))
```

In addition, the [Pagoda2 pacakge](https://github.com/kharchenkolab/pagoda2) has to be installed to proces. Please check on the corresponding github page the different dependencies needed to install it. 

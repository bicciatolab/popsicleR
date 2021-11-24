# _PoPsicleR_ 
### A flexible R Package for the pre-processing and quality control analysis of single cell RNA-seq data
__________________________________________________________________
#### Summary: 
Here we present PoPsicleR, an R package for the preprocessing and QC of single cell RNA-seq data that allows users to interactively optimize thresholds and parameters of the commonly used QC metrics. PoPsicleR builds on the preprocessing functions of the Seurat package (v3.5) for filtering, normalization, regression, clustering, annotation, and data visualization and integrates the Python module of Scrublet for doublet detection. As compared to other preprocessing and QC pipelines, PoPsicle offers specific functions to optimize filtering thresholds, design user-specific preprocessing setups, and guide inexperienced command line users through all the various steps of the preprocessing workflow. PoPsicle main assets are its efficacy in testing complex combinations of filtering thresholds in a simple manner and its lightweight computational requirements that enable scientists to process large number of cells obtained from different technologies (as 10X and Smart-seq) and species (e.g., human and mouse).

#### Availability and implementation:
PoPsicleR is written in R language and is released under a GPL License. It can be downloaded from Github ([@bicciatolab/PoPsicleR](https://github.com/bicciatolab/PoPsicleR))
#### Contact: 
silvio.bicciato@unimore.it


## Table of contents

- `PoPsicleR` [installation](https://github.com/bicciatolab/PoPsicleR/blob/main/docs/PoPsicleR_installation_guide.md#popsicler-installation-in-r) in R
- `PoPsicleR` [installation](https://github.com/bicciatolab/PoPsicleR/blob/main/docs/PoPsicleR_installation_guide.md#popsicler-installation-through-conda) through `conda`
- [Quick Start guide](https://github.com/bicciatolab/PoPsicleR/docs/Quick_Start_guide.md)

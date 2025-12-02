# popsicleR

<img src="https://github.com/bicciatolab/popsicleR/blob/main/man/figures/Logo_popsicleR.png" width=20% height= 20% alt="Logo_popsicleR"  align= "right">

`popsicleR` is a R package that combines methods implemented in widely used pipelines to interactively perform all major pre-processing and QC steps of scRNA-seq data analysis. The package is composed of seven functions capable of performing exploration of quality-control metrics, filtering of low-quality cells, data normalization, removal of technical and biological biases, and some basic analysis as detection of differentially expressed genes, cell clustering and cell annotation. During each step of the analysis, `popsicleR` interactively guides the user with colored text messages and saves in dedicated folders a variety of plots to investigate several QC metrics and assess the impact of filtering and regression parameters on the identification and classification of cell populations.

Key features of `popsicleR` include:
1. Use as input of files from either the Cell Ranger pipeline of 10X Genomics or a feature-barcode matrix of raw counts generated from any microfluidic-, microwell plates-, or droplet-based scRNA-seq technology
2. Output of graphs and colored text messages to interactively guide users along each step of the analysis
3. Inclusion of common single-cell visualisations (as density, scatter, and violin plots and low-dimensionality embeddings) to investigate QC metrics and pre-processing parameters
4. Export of visualisations as PDF images for presentation or publication use

<p align="center">
<img src="https://github.com/bicciatolab/popsicleR/blob/main/docs/popsicleR_workflow.png" width="40%" alt="popsicleR-workflow">
</p>

#### Contact:

silvio.bicciato@unimore.it; mattia.forcato@unimore.it

#### Citation:

F. Grandi, J. Caroli, O. Romano, M. Marchionni, M. Forcato, S. Bicciato, popsicleR: a R Package for pre-processing and quality control analysis of single cell RNA-seq data, _Journal of Molecular Biology_ (2022),  doi: [10.1016/j.jmb.2022.167560](https://doi.org/10.1016/j.jmb.2022.167560)

# Table of Contents

- [System requirements](https://github.com/bicciatolab/popsicleR#System-requirements)
- [Installation](https://github.com/bicciatolab/popsicleR#installation)
- [Tutorial on example data](https://raw.githack.com/bicciatolab/popsicleR/main/docs/popsicleR_tutorial.html)

## System requirements

* R version: >= 4.0.0
* Dependencies: *ape*, *celldex*, *clustree*, *corrplot*, *crayon*, *dplyr*, *future*, *ggExtra*, *ggplot2*, *ggplotify*, *gtools*, *grid*, *gridExtra*, *limma*, *magrittr*, *patchwork*, *pheatmap*, *neldermead*, *RANN*, *RColorBrewer*, *reticulate*, *R.utils*, *scDblFinder*, *scMCA*, *session*, *shinythemes*, *umap*, *Seurat*, and *SingleR*.

## Installation

In order to avoid conflicts between package dependencies we provide here a comprehensive guide to install `popsicleR` through Anaconda platform.

If Anaconda is already installed, a virtual environment for `popsicleR` can be set either manually, installing all packages one by one, or automatically adding only selected packages once created the environment.

The following comands allow setting the `popsicleR` environment.

#### Create a `popsicleR` environment and install all anaconda available packages automatically

To create the `popsicleR` environment on a Linux machine, open the terminal and run:

```bash
conda create -n popsicleR -c conda-forge r-base r-umap r-neldermead r-rann r-rcolorbrewer r-ggextra r-ggplotify r-crayon r-patchwork r-magrittr r-gridextra r-dplyr r-ggplot2 r-devtools r-r.utils r-future r-reticulate r-pheatmap r-shinythemes r-rcurl r-corrplot r-locfit r-clustree r-ape r-sessioninfo r-seurat bioconda::bioconductor-singler bioconda::bioconductor-limma bioconda::bioconductor-celldex bioconda::bioconductor-scdblfinder r::r-session r-matrixstats=1.1.0 r-igraph=1.5.0
```

#### Install `popsicleR` environment through a .yml file

`conda popsicleR` environment can also be extracted from a [popsicleR.yml](https://github.com/bicciatolab/popsicleR/blob/main/docs/popsicleR.yml) file. In this case, all anaconda required packages will be automatically installed. After downloading the [popsicleR.yml](https://github.com/bicciatolab/popsicleR/blob/main/docs/popsicleR.yml) file in the working directory (or on a specific file_path), run:

```bash
conda env create -n popsicleR -f popsicleR.yml
```

#### Install popsicleR package from Github

Once created the environment, access it through the command:

```bash
conda activate popsicleR
```

Before proceeding with the installation of popsicleR, scMCA package must be installed from Github repository.
Thus, open `R` and use the following scripts to install `scMCA` package from Github:

```r
install.packages("stringi")
devtools::install_github("ggjlab/scMCA")
```

Then, **popsicleR** package can be, similarly, installed from Github:

```r
devtools::install_github("bicciatolab/popsicleR")
```

In case of any issue with installation of `popsicleR` via `install_github`, it is possible to download the package.tar.gz from [here](https://github.com/bicciatolab/popsicleR/popsicleR_0.3.0.tar.gz) using the bash command:

```bash
 wget https://github.com/bicciatolab/popsicleR/popsicleR_0.3.0.tar.gz
```

Lately, extract the main directory (if necessary, rename the package folder from "popsicleR-main" to "popsicleR") and install `popsicleR` from the local repository with the following script: 

```r
install.packages("/path/to/package_directory", repos = NULL, type="source")
```

#### Backward compatibility with Seurat version 4

Recently, **popsicleR** was updated to work with Seurat version 5.X, to ensure compatibility with the newest Seurat releases. If the old version of Seurat 4 is desired, install the conda environment using a different .yml file. After downloading the [popsicleR_Seurat4.yml](https://github.com/bicciatolab/popsicleR/blob/main/docs/popsicleR_Seurat4.yml) file in the working directory (or on a specific file_path), run:

```bash
conda env create -n popsicleR -f popsicleR_Seurat4.yml
```

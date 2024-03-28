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
conda create -n popsicleR -c conda-forge r-base=4.0.3 r-umap=0.2.7.0 r-neldermead=1.0_11 r-rann=2.6.1 r-rcolorbrewer=1.1_2 r-ggextra=0.9 r-ggplotify=0.1.0 r-crayon=1.4.0 r-patchwork=1.1.1 r-magrittr=1.5 r-gridextra=2.3 r-dplyr=1.0.4 r-ggplot2=3.3.3 r-devtools=2.3.2 r-r.utils=2.10.1 r-future=1.21.0 r-reticulate=1.18 r-pheatmap=1.0.12 r-shinythemes=1.2.0 r-rcurl=1.98_1.2 r-corrplot=0.92 r-locfit=1.5_9.4 r-clustree=0.4.4 r-ape=5.6 r-seuratobject=4.0.4 r-sessioninfo=1.1.1 r-seurat
```

#### Install environment packages

From command line, use the following `conda`  commands to install packages from other channels:

```bash
conda install -n popsicleR -c r r-magrittr
conda install -n popsicleR -c bioconda bioconductor-limma=3.46.0
```

Since not all required packages are provided in anaconda.org, some packages must be intalled directly from the R console, as described in [install packages in R](https://github.com/bicciatolab/popsicleR#install-packages-in-r).

#### Install `popsicleR` environment through a .yml file

`conda popsicleR` environment can also be extracted from a [popsicleR.yml](https://github.com/bicciatolab/popsicleR/blob/main/docs/popsicleR.yml) file. In this case, all anaconda required packages will be automatically installed. After downloading the [popsicleR.yml](https://github.com/bicciatolab/popsicleR/blob/main/docs/popsicleR.yml) file in the working directory (or on a specific file_path), run:

```bash
conda env create -n popsicleR -f popsicleR.yml
```

Since not all required packages are provided in anaconda.org, some packages must be intalled directly from the R console, as described in [install packages in R](https://github.com/bicciatolab/popsicleR#install-packages-in-r).

#### Install packages in R
Once created the environment, access it through the command:

```bash
conda activate popsicleR
```
Before installing `popsicleR`, users must run the following codes to install packages required as dependencies. 
During this stage, we recommend skipping updates and the installation of `Rtools` when prompted by R. 

Thus, open `R` and install the dependencies `SingleR`, `celldex`, `scDblFinder` from Bioconductor and `scMCA`  from Github using:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scDblFinder")
```

```r
devtools::install_github("ggjlab/scMCA")
```

Before the last step, install `session` and `dbplyr` packages directly from CRAN archive using devtools.

```r
devtools::install_version("session", version = "1.0.3", repos = "https://cran.r-project.org/")
devtools::install_version("dbplyr", version = "2.3.4", repos = "https://cran.r-project.org/")
```

Finally, use the following scripts to install `popsicleR` from Github :

```r
devtools::install_github("bicciatolab/popsicleR")
```
In case of any issue with installation of `popsicleR` via `install_github`, it is possible to download the package.tar.gz using the bash command:

```bash
 wget https://github.com/bicciatolab/popsicleR/archive/main.tar.gz
```

Lately, extract the main directory (if necessary, rename the package folder from "popsicleR-main" to "popsicleR") and install `popsicleR` from the local repository with the following script:

```r
install.packages("/path/to/package_directory", repos = NULL, type="source")
```
# _PopsicleR_ 
### A flexible R Package for the pre-processing and quality control analysis of single cell RNA-seq data
__________________________________________________________________
#### Summary: 
Here we present `popsicleR`, an R package for the preprocessing and QC of single cell RNA-seq data that allows users to interactively optimize thresholds and parameters of the commonly used QC metrics. `popsicleR` builds on the preprocessing functions of the Seurat package for filtering, normalization, regression, clustering, annotation, and data visualization and integrates the Python module of Scrublet for doublet detection. As compared to other preprocessing and QC pipelines, `popsicleR` offers specific functions to optimize filtering thresholds, design user-specific preprocessing setups, and guide inexperienced command line users through all the various steps of the preprocessing workflow. `popsicleR` main assets are its efficacy in testing complex combinations of filtering thresholds in a simple manner and its lightweight computational requirements that enable scientists to process large number of cells obtained from different technologies (as 10X and Smart-seq) and species (e.g., human and mouse).

#### Availability and implementation:
popsicleR is written in R language and is released under a GPL License. It can be downloaded from Github ([@bicciatolab/popsicleR](https://github.com/bicciatolab/popsicleR))

#### Contact: 
silvio.bicciato@unimore.it


## Table of contents

- `popsicleR` [installation in R](https://github.com/bicciatolab/popsicleR#popsicler-installation-in-r) 
- `popsicleR` [installation through `conda`](https://github.com/bicciatolab/popsicleR#popsicler-installation-through-conda) 
- [Quick Start guide](https://github.com/bicciatolab/popsicleR/docs/Quick_Start_guide.md)

## popsicleR installation in R

Open your R and copy and paste on your console these lines of code to install all the required packages: 

```
install.packages("lattice")
install.packages("session")

## install packages specific versions from source:
   
packagesurl <- c("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-2.tar.gz","https://cran.r-project.org/src/contrib/Archive/optimbase/optimbase_1.0-9.tar.gz","https://cran.r-project.org/src/contrib/Archive/optimsimplex/optimsimplex_1.0-7.tar.gz","https://cran.r-project.org/src/contrib/Archive/neldermead/neldermead_1.0-11.tar.gz")

install.packages(packagesurl, repos=NULL, type="source")

## Seurat

install.packages("Seurat")

## SingleR

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleR")

## devtools

install.packages("devtools")
library("devtools")

## celldex and limma: 

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("celldex")
BiocManager::install("limma")

## pheatmap, scMCA, popsicleR:

install_github("cran/pheatmap") 
install_github("ggjlab/scMCA") 
install_github("bicciatolab/popsicleR")
```
or download the package.tar.gz file and then install it from local repository through the command:
  
```
install.packages("/path/to/package_directory", repos = NULL, type="source")
```

## popsicleR installation through conda

If you have already installed Anaconda, a virtual environment for popsicleR can be set either manually, installing all packages one by one, or automatically adding only a few packages once created the environment.
Another way to set an environment is to extract it from a file.yml, we provide popsicle.yml file; in this case all anaconda required packages will be automatically installed.

Below we report all the instruction to set your popsicleR environment.

#### Create a popsicleR environment and install all anaconda available packages automatically
 
On a Linux machine, open your terminal and run:

```bash
	$conda create -n popsicleR -c conda-forge r-base=4.0.3 r-umap=0.2.7.0 r-neldermead=1.0_11 r-rann=2.6.1 r-rcolorbrewer=1.1_2 r-ggextra=0.9 r-ggplotify=0.1.0 r-crayon=1.4.0 r-patchwork=1.1.1 r-magrittr=1.5 r-gridextra=2.3 r-dplyr=1.0.4 r-ggplot2=3.3.3 r-devtools=2.3.2 r-r.utils=2.10.1 r-future=1.21.0 r-reticulate=1.18 r-pheatmap=1.0.12 r-shinythemes=1.2.0 r-rcurl=1.98_1.2 r-seuratobject=4.0.0 r-sessioninfo=1.1.1 r-seurat
```
#### install environment packages

install some other packages from other channels through conda command: 

```bash
	$conda install -n popsicleR -c r r-magrittr
	$conda install -n popsicleR -c bioconda bioconductor-limma=3.46.0
	$conda install -n popsicleR -c bioconda bioconductor-singler
	$conda install -n popsicleR -c bioconda bioconductor-celldex
```
Then scroll down and go to the section "install packages in R"

#### install popsicleR environment through a .yml file

Be secure that "popsicleR.yml" is present in your working directory then run: 

```bash
	$conda env create -n popsicleR -f popsicleR.yml
```

Then scroll down and go to the section "install packages in R"

#### install packages in R

Pay attention to the released packages on anaconda.org; if needed packages are provided in anaconda repository we suggests to install it through conda command in order to avoid conflicts between packages.

Step 2: access the environment 

Once created the environment and installed all the available packages like explained before, you can access the popsicleR virtual environment with the command: 

```bash
	$conda activate popsicleR
```
The last dependencies, not provided at the moment by anaconda.org, can be installed directly from R.

```r
library(devtools)
install_github("ggjlab/scMCA") 

# N.B. if "ggplotify" is not v0.1.0 update the package, and then run:

devtools::install_github("bicciatolab/popsicleR")
```
or download the package.tar.gz file and then install it from local  repository through the command:
 
```
install.packages("/path/to/package_directory", repos = NULL, type="source")
```
# _PoPsicleR_ 
### A flexible R Package for the pre-processing and quality control analysis of single cell RNA-seq data
__________________________________________________________________
#### Summary: 
Here we present `PoPsicleR`, an R package for the preprocessing and QC of single cell RNA-seq data that allows users to interactively optimize thresholds and parameters of the commonly used QC metrics. `PoPsicleR` builds on the preprocessing functions of the Seurat package (v3.5) for filtering, normalization, regression, clustering, annotation, and data visualization and integrates the Python module of Scrublet for doublet detection. As compared to other preprocessing and QC pipelines, `PoPsicleR` offers specific functions to optimize filtering thresholds, design user-specific preprocessing setups, and guide inexperienced command line users through all the various steps of the preprocessing workflow. `PoPsicleR` main assets are its efficacy in testing complex combinations of filtering thresholds in a simple manner and its lightweight computational requirements that enable scientists to process large number of cells obtained from different technologies (as 10X and Smart-seq) and species (e.g., human and mouse).

#### Availability and implementation:
PoPsicleR is written in R language and is released under a GPL License. It can be downloaded from Github ([@bicciatolab/PoPsicleR](https://github.com/bicciatolab/PoPsicleR))
#### Contact: 
silvio.bicciato@unimore.it


## Table of contents

- `PoPsicleR` [installation](https://github.com/bicciatolab/PoPsicleR#table-of-contents) in R
- `PoPsicleR` [installation](https://github.com/bicciatolab/PoPsicleR/blob/main/docs/PoPsicleR_installation_guide_through_conda.md) through `conda`
- [Quick Start guide](https://github.com/bicciatolab/PoPsicleR/docs/Quick_Start_guide.md)

## PoPsicleR installation in R

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

## pheatmap, scMCA, popsicler:

install_github("cran/pheatmap") 
install_github("ggjlab/scMCA") 
install_github("bicciatolab/PoPsicleR")
```
or download the package.tar.gz file and then install it from local repository through the command:
  
```
install.packages("/path/to/package_directory", repos = NULL, type="source")
```

## PoPsicleR installation through conda

If you have already installed Anaconda, a virtual environment for PoPsicleR can be set either manually, installing all packages one by one, or automatically adding only a few packages once created the environment.
Another way to set an environment is to extract it from a file.yml, we provide PoPsicle.yml file; in this case all anaconda required packages will be automatically installed.

Below we report all the instruction to set your PoPsicleR environment.


##### Create a PoPsicle environment and install all anaconda available packages manually 

On a Linux machine, open your terminal and run:

Step1: create a conda environment

	$conda create -n PoPsicleR

Step 2: access the environment 

	$conda activate PoPsicleR

Step 3: manual installation

	$conda install -c conda-forge r-base     	#R 4.0.3
	$conda install -c conda-forge r-umap		#umap 0.2.7.0
	$conda install -c conda-forge r-neldermead	#neldermead 1.0_11
	$conda install -c conda-forge r-rann		#RANN 2.6.1
	$conda install -c conda-forge r-rcolorbrewer	#Rcolorbrewer 1.1_2
	$conda install -c conda-forge r-ggextra		#ggExtra 0.9
	$conda install -c conda-forge r-ggplotify	#ggplotify 0.0.5
	$conda install -c conda-forge r-crayon		#crayon 1.4.0
	$conda install -c conda-forge r-patchwork	#patchwork 1.1.1
	$conda install -c bioconda bioconductor-limma	#limma 3.46.0
	$conda install -c r r-magrittr			#magrittr 1.5
	$conda install -c conda-forge r-gridextra	#gridExtra 2.3
	$conda install -c conda-forge r-dplyr		#dplyr 1.0.4
	$conda install -c conda-forge r-ggplot2		#ggplot2 3.3.3
	$conda install -c conda-forge r-devtools	#devtools 2.3.2
	$conda install -c conda-forge r-r.utils		#r.utils 2.10.1
	$conda install -c conda-forge r-future		#future 1.21.0
	$conda install -c conda-forge r-reticulate	#reticulate 1.18
	$conda install -c conda-forge r-pheatmap	#pheatmap 1.0.12
	$conda install -c conda-forge r-shinythemes	#shinythemes 1.2.0
	$conda install -c conda-forge r-rcurl		#RCurl 1.98_1.2
	$conda install -c conda-forge r-seuratobject	#SeuratObject 4.0.0

Then scroll down and go to the section "install packages in R"

##### Install PoPsicle environment automatically
 
On a Linux machine, open your terminal and run:

	$conda create -n PoPsicleR -c conda-forge r-base=4.0.3 r-umap=0.2.7.0 r-neldermead=1.0_11 r-rann=2.6.1 r-rcolorbrewer=1.1_2 r-ggextra=0.9 r-ggplotify=0.1.0 r-crayon=1.4.0 r-patchwork=1.1.1 r-magrittr=1.5 r-gridextra=2.3 r-dplyr=1.0.4 r-ggplot2=3.3.3 r-devtools=2.3.2 r-r.utils=2.10.1 r-future=1.21.0 r-reticulate=1.18 r-pheatmap=1.0.12 r-shinythemes=1.2.0 r-rcurl=1.98_1.2 r-seuratobject=4.0.0 r-sessioninfo=1.1.1

##### install environment packages

install some other packages from other channels through conda command: 

	$conda install -n PoPsicleR -c r r-magrittr
	$conda install -n PoPsicleR -c bioconda bioconductor-limma=3.46.0		#limma 3.46.0

Then scroll down and go to the section "install packages in R"

##### install PoPsicle environment through a .yml file

Be secure that "PoPsicle.yml" is present in your working directory then run: 

	$conda env create -n PoPsicleR -f PoPsicleR.yml

Then scroll down and go to the section "install packages in R"

#####	install packages in R (Required for each installation method)   

Pay attention to the released packages on anaconda.org; if needed packages are provided in anaconda repository we suggests to install it through conda command in order to avoid conflicts between packages.

Once created the environment and installed all the available packages like explained before, you can access the PoPsicleR virtual environment with the command: 

	$conda activate PoPsicleR

The last dependencies, not provided at the moment by anaconda.org, can be installed directly from R.
 
```
install.packages('Seurat')

# SingleR

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleR")

# celldex,  scMCA,  popsicler:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("celldex")

library(devtools)
install_github("ggjlab/scMCA") 

# N.B. if "ggplotify" is not v0.1.0 update the package, and then run:

devtools::install_github("bicciatolab/PoPsicleR")
```
or download the package.tar.gz file and then install it from local  repository through the command:
 
```
install.packages("/path/to/package_directory", repos = NULL, type="source")
```


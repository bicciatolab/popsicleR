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

__________________________________________________________________
## popsicleR installation in R
 
On a Windows or Linux machine, open your R/Rstudio and install all the required packages using these lines of code.

```r
install.packages("devtools")
library("devtools")
install_github("cran/pheatmap") 
install_github("ggjlab/scMCA") 
install_github("bicciatolab/popsicleR")
```
`install_github` by default try to install all dependencies; but sometimes can be needed to manually install packages from specific repository. 

In this cases, as a first step, check if you have already install all/some or none of these specific packages, through:

```r
requiredPackages <- c("SingleR","limma","neldermead","celldex","Matrix","optimbase","optimsimplex") 

newPackages <- requiredPackages[!(requiredPackages %in% installed.packages()[,"Package"])]
```
`newPackages` object return a list of not already installed packages. 

To install these new packages, run:

```r
specificversionDependencies <- c("Matrix","optimbase","optimsimplex","neldermead")

if (length(newPackages[newPackages %in% specificversionDependencies])>0) 
	{
	to_install <- newPackages[newPackages %in% specificversionDependencies]
	reorder <- match(specificversionDependencies, to_install)
	to_install <- to_install[reorder]
	to_install <- to_install[!is.na(to_install)]
	packagesurl <- c("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-2.tar.gz","https://cran.r-project.org/src/contrib/Archive/optimbaseoptimbase_1.0-9.tar.gz","https://cran.r-project.org/src/contrib/Archive/optimsimplex/optimsimplex_1.0-7.tar.gz","https://cran.r-project.org/src/contrib/Archive/neldermead/neldermead_1.0-11.tar.gz")
	for (i in 1:length(to_install))	
		{
		source_repo <- packagesurl[grep(to_install[i], packagesurl)]
		install.packages(source_repo, repos=NULL, type="source")
		} 


BiocManagerDependencies <- c("SingleR","limma", "celldex")

if (length(newPackages[newPackages %in% BiocManagerDependencies])>0) 
	{
		to_install <- newPackages[newPackages %in% BiocManagerDependencies]
		reorder <- match(BiocManagerDependencies, to_install)
		to_install <- to_install[reorder]
		to_install <- to_install[!is.na(to_install)]
		if (!requireNamespace("BiocManager", quietly = TRUE))
    		install.packages("BiocManager")
    	BiocManager::install(to_install)
	} 
```

Now you should be able to install `popsicleR` via: 

```r
install_github("bicciatolab/popsicleR")
```
If you had problem with this installation procedure download the package.tar.gz and go to the  [popsicleR manual installation](https://github.com/bicciatolab/popsicleR/tree/main/docs) page: there you will find all the required packages dependencies managed and specified. Finally, install `popsicleR` from a local repository through the command:
```r
install.packages("/path/to/package_directory", repos = NULL, type="source")
```

## popsicleR installation through conda

If you have already installed Anaconda, a virtual environment for popsicleR can be set either manually, installing all packages one by one, or automatically adding only a few packages once created the environment.

Below we report all the instruction to set your `popsicleR` environment.

#### Create a popsicleR environment and install all anaconda available packages automatically
 
`popsicleR` environment creation
  
On a Linux machine, open your terminal and run:

```bash
	conda create -n popsicleR -c conda-forge r-base=4.0.3 r-umap=0.2.7.0 r-neldermead=1.0_11 r-rann=2.6.1 r-rcolorbrewer=1.1_2 r-ggextra=0.9 r-ggplotify=0.1.0 r-crayon=1.4.0 r-patchwork=1.1.1 r-magrittr=1.5 r-gridextra=2.3 r-dplyr=1.0.4 r-ggplot2=3.3.3 r-devtools=2.3.2 r-r.utils=2.10.1 r-future=1.21.0 r-reticulate=1.18 r-pheatmap=1.0.12 r-shinythemes=1.2.0 r-rcurl=1.98_1.2 r-seuratobject=4.0.4 r-sessioninfo=1.1.1 r-seurat
```

#### Install environment packages

From command line, install packages from other channels using ` conda`  commands: 

```bash
	conda install -n popsicleR -c r r-magrittr
	conda install -n popsicleR -c bioconda bioconductor-limma=3.46.0 
```

Not all the required packages are provided in anaconda.org and the last packages must be intalled directly from R consolle: 

Please go to the [install packages in R](https://github.com/bicciatolab/popsicleR#install-packages-in-r) section to complete `popsicleR` installation.

#### Install popsicleR environment through a .yml file

Another way to set a working `conda` popsicleR environment is to  extract it from a .yml file. 

In this case, all anaconda required packages will be automatically installed.

To easily create your `popsicleR` environment, we provide a [popsicleR.yml](https://github.com/bicciatolab/popsicleR/blob/main/docs/popsicleR.yml) file in the [docs](https://github.com/bicciatolab/popsicleR/tree/main/docs) folder. 

Thus, download [popsicleR.yml](https://github.com/bicciatolab/popsicleR/blob/main/docs/popsicleR.yml) file, be secure that it is present in your working directory (or specify file_path) and finally, run: 

```bash
	conda env create -n popsicleR -f popsicleR.yml
```

Not all the required packages are provided in anaconda.org and the last packages must be intalled directly from R consolle: 

Please go to the [install packages in R](https://github.com/bicciatolab/popsicleR#install-packages-in-r) section to complete `popsicleR` installation.

#### Install packages in R

The last dependencies, at the moment not provided by anaconda.org, can be installed directly from R.

Once created the environment, access it through the command: 

```bash
	conda activate popsicleR
```

install `SingleR`, `celldex` and `scMCA` packages using:


```r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleR")
BiocManager::install("celldex")

devtools::install_github("ggjlab/scMCA") 
```

and finally, install `popsicleR` via:

```r
devtools::install_github("bicciatolab/popsicleR")
```

If you have with this command download the package.tar.gz file and then install it from local repository through the command:
 
```r
install.packages("/path/to/package_directory", repos = NULL, type="source")
```

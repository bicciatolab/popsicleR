# PopsicleR 

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

# Table of Contents

- [System requirements](https://github.com/bicciatolab/popsicleR#System-requirements)
- [Installation in R](https://github.com/bicciatolab/popsicleR#installation-in-r)
- [Installation through `conda`](https://github.com/bicciatolab/popsicleR#installation-through-conda) 
- [Tutorial on example data](http://htmlpreview.github.io/?https://github.com/bicciatolab/popsicleR/main/docs/popsicleR_tutorial.html)

## System requirements

* R version: >= 4.0.0
* Dependencies: *ape*, *celldex*, *clustree*, *corrplot*, *crayon*, *dplyr*, *future*, *ggExtra*, *ggplot2*, *ggplotify*, *gtools*, *grid*, *gridExtra*, *limma*, *magrittr*, *patchwork*, *pheatmap*, *neldermead*, *RANN*, *RColorBrewer*, *reticulate*, *R.utils*, *scMCA*, *session*, *shinythemes*, *umap*, *Seurat*, and *SingleR*.

## Installation in R

Before installing `popsicleR`, users can run the following codes to install packages from CRAN, CRAN archive, Bioconductor, and Github required as dependencies:

```r
CRANdep <- c("Seurat","reticulate","R.utils","dplyr","ggplot2","clustree","ape","gtools",
"future","grid","gridExtra","magrittr","limma","patchwork",
"crayon","ggExtra","RColorBrewer","ggplotify","RANN","umap",
"celldex","curl","httr","lattice","session","shinythemes","usethis","rcmdcheck",
"roxygen2","rversions","devtools","pheatmap","BiocManager")
newPackages <- CRANdep[!(CRANdep %in% installed.packages()[,"Package"])]
if(length(newPackages)){install.packages(newPackages)}

CRANarcdep <- c("Matrix","optimbase","optimsimplex","neldermead")
newPackages <- CRANarcdep[!(CRANarcdep %in% installed.packages()[,"Package"])]
if(length(newPackages)){
	packagesurl <- c("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-2.tar.gz",
	"https://cran.r-project.org/src/contrib/Archive/optimbase/optimbase_1.0-9.tar.gz",
	"https://cran.r-project.org/src/contrib/Archive/optimsimplex/optimsimplex_1.0-7.tar.gz",
	"https://cran.r-project.org/src/contrib/Archive/neldermead/neldermead_1.0-11.tar.gz")
	for (i in 1:length(newPackages)){
		source_repo <- packagesurl[grep(newPackages[i], packagesurl)]
		install.packages(source_repo, repos=NULL, type="source")
	} 
}

BioCdep <- c("SingleR","limma","BiocFileCache","AnnotationHub","ExperimentHub", "celldex")
newPackages <- BioCdep[!(BioCdep %in% installed.packages()[,"Package"])]
if(length(newPackages)){BiocManager::install(newPackages)}

if(!"scMCA"%in% installed.packages()[,"Package"]){devtools::install_github("ggjlab/scMCA")}

```
 
Once installed all dependencies, you can install `popsicleR` with the following script:

```r
devtools::install_github("bicciatolab/popsicleR")
```

In case of any issue with installation of `popsicleR` via `install_github`, it is possible to download the package.tar.gz and install `popsicleR` from the local repository with the following script:

```r
install.packages("/path/to/package_directory", repos = NULL, type="source")
```

## Installation through conda

If Anaconda is already installed, a virtual environment for `popsicleR` can be set either manually, installing all packages one by one, or automatically adding only selected packages once created the environment.

The following comands allow setting the `popsicleR` environment.

#### Create a `popsicleR` environment and install all anaconda available packages automatically
  
To create the `popsicleR` environment on a Linux machine, open the terminal and run:

```bash
conda create -n popsicleR -c conda-forge r-base=4.0.3 r-umap=0.2.7.0 r-neldermead=1.0_11 r-rann=2.6.1 r-rcolorbrewer=1.1_2 r-ggextra=0.9 r-ggplotify=0.1.0 r-crayon=1.4.0 r-patchwork=1.1.1 r-magrittr=1.5 r-gridextra=2.3 r-dplyr=1.0.4 r-ggplot2=3.3.3 r-devtools=2.3.2 r-r.utils=2.10.1 r-future=1.21.0 r-reticulate=1.18 r-pheatmap=1.0.12 r-shinythemes=1.2.0 r-rcurl=1.98_1.2 r-seuratobject=4.0.4 r-sessioninfo=1.1.1 r-seurat
```

#### Install environment packages

From command line, use the following `conda`  commands to install packages from other channels: 

```bash
conda install -n popsicleR -c r r-magrittr
conda install -n popsicleR -c bioconda bioconductor-limma=3.46.0 
```

Since not all required packages are provided in anaconda.org, some packages must be intalled directly from the R console, as described in [install packages in R](https://github.com/bicciatolab/popsicleR#install-packages-in-r).

#### Install `popsicleR` environment through a .yml file

`conda popsicleR` environment can also be extracted from a [popsicleR.yml file](https://github.com/bicciatolab/popsicleR/blob/main/docs/popsicleR.yml). In this case, all anaconda required packages will be automatically installed. After downloading the [popsicleR.yml file](https://github.com/bicciatolab/popsicleR/blob/main/docs/popsicleR.yml) in the working directory (or on a specific file_path), run: 

```bash
conda env create -n popsicleR -f popsicleR.yml
```

Since not all required packages are provided in anaconda.org, some packages must be intalled directly from the R console, as described in [install packages in R](https://github.com/bicciatolab/popsicleR#install-packages-in-r).

#### Install packages in R
Once created the environment, access it through the command: 

```bash
conda activate popsicleR
```

and install `SingleR`, `celldex` and `scMCA` packages using:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleR")
BiocManager::install("celldex")

devtools::install_github("ggjlab/scMCA") 
```

Finally, use the following scripts to install `popsicleR` from Github :

```r
devtools::install_github("bicciatolab/popsicleR")
```

or from a local repository:
 
```r
install.packages("/path/to/package_directory", repos = NULL, type="source")
```
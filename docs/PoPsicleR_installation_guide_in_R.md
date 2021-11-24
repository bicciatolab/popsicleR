# PoPsicleR installation in R

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
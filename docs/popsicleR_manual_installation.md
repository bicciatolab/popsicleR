If you had problem with the installation procedure, download the package.tar.gz and install all the required packages as described below. 

As a first step, open your R/Rstudio and check if you have already install all/some or none of these needed packages.

```r
requiredPackages <- c("Seurat","reticulate","R.utils","dplyr","ggplot2","clustree","ape","gtools","SingleR","future","session","grid","gridExtra","magrittr","limma","patchwork","crayon","ggExtra","RColorBrewer","ggplotify","RANN","neldermead","umap","celldex","scMCA","devtools","Matrix","optimbase","optimsimplex","curl","httr","BiocFileCache","AnnotationHub","ExperimentHub","lattice","session","pheatmap")

newPackages <- requiredPackages[!(requiredPackages %in% installed.packages()[,"Package"])]
```
Then run:

```r
CRANDependencies <- c("Seurat","reticulate","R.utils","dplyr","ggplot2","clustree","ape","gtools","future","grid","gridExtra","magrittr","limma","patchwork","crayon","ggExtra","RColorBrewer","ggplotify","RANN","umap","celldex","curl","httr","lattice","session","usethis","rcmdcheck","roxygen2","rversions","devtools")

if (length(newPackages[newPackages %in% CRANDependencies])>0) {
	to_install <- newPackages[newPackages %in% CRANDependencies]
	reorder <- match(CRANDependencies, to_install)
	to_install <- to_install[reorder]
	to_install <- to_install[!is.na(to_install)]
	install.packages(to_install)
} 

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

BiocManagerDependencies <- c("SingleR","limma","BiocFileCache","AnnotationHub","ExperimentHub", "celldex")

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

githubDependencies <- c("pheatmap","scMCA")

if (length(newPackages[newPackages %in% githubDependencies])>0) 
    {
		to_install <- newPackages[newPackages %in% githubDependencies]
		reorder <- match(githubDependencies, to_install)
		to_install <- to_install[reorder]
		to_install <- to_install[!is.na(to_install)]
		git_sources <- c("cran/pheatmap","ggjlab/scMCA")
		for (i in 1:length(to_install)){
		source_repo <- git_sources[grep(to_install[i], git_sources)]
		devtools::install_github(source_repo)
	} 
```

Taking advantage of these if statement, you should be able to manage all the packages coming from different repositories and to set a working environment for `popsicleR`.

Lastly, install `popsicleR` from a local repository through the command:

```r
install.packages("/path/to/package_directory", repos = NULL, type="source")
```
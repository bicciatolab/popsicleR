\name{PoPsicle}
\alias{popsicle}
\alias{popsicle}
\docType{package}
\title{
\packageTitle{popsicle}
}
\description{

\packageDescription{popsicle}


}
\details{

PoPsicle builds on the R package Seurat (v3.5) to perform most of its task, and integrates an R implementation of the Scrublet python module for the analysis of doublets presence. Furthermore, PoPsicle grants the possibility to investigate and analyze both human and mouse data, using organism specific integrated annotation packages.

\figure{popsicle_logo.png}{options: width=200 alt "Figure: popsicle_logo.png" align= right}


\bold{Functions Index:}
\itemize{
\item \code{\link{PrePlots}}
\item \code{\link{FilterPlots}}
\item \code{\link{CalculateDoublets}}
\item \code{\link{Normalize}}
\item \code{\link{CCScore}}
\item \code{\link{ApplyRegression}}
\item \code{\link{CalculateCluster}}
\item \code{\link{MakeAnnotation}}
}

}
\author{
\packageAuthor{popsicle}

\bold{Maintainers}:
\itemize{
\item Jimmy Caroli \email{jimmy.caroli@sund.ku.dk}
\item Francesco Grandi \email{francesco.grandi1993@unimore.it}
}
}

\references{
Stuart T, Butler A, Hoffman P, Hafemeister C, Papalexi E, III WMM, Hao Y, Stoeckius M, Smibert P, Satija R (2019). “Comprehensive Integration of Single-Cell Data.” \emph{Cell}, \bold{177}, --1888-1902--. \url{https://doi.org/10.1016/j.cell.2019.05.031.}

Samuel L. Wolock, Romain Lopez, Allon M. Klein, "Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data" (2019). \emph{Cell System}, Volume \bold{8},--281-291--.e9,ISSN 2405-4712,\url{https://doi.org/10.1016/j.cels.2018.11.005.}

Aran D, Looney AP, Liu L, Wu E, Fong V, Hsu A, Chak S, Naikawadi RP, Wolters PJ, Abate AR, Butte AJ, Bhattacharya M (2019). “Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.” \emph{Nat. Immunol.}, \bold{20}, --163-172--. \url{doi: 10.1038/s41590-018-0276-y.}

Sun H, Zhou Y, Fei L, Chen H, Guo G. "scMCA: A Tool to Define Mouse Cell Types Based on Single-Cell Digital Expression." (2019) \emph{Methods Mol Biol.}, \bold{1935}, --91-96--. \url{doi: 10.1007/978-1-4939-9057-3_6. PMID: 30758821.}

}
\keyword{ popsicle, scRNA-seq, Seurat, SingleR, scMCA }
\seealso{
Useful links:
\itemize{
    \item \url{https://github.com/bicciatolab/PoPsicle}

    \item Report bugs at \url{https://github.com/bicciatolab/PoPsicle/issues}
}
}
\examples{
PrePlots("sample_name", input_data= "/path/to/cellranger/output/folder", genelist=c("PTPRC","CD34","PROM1","CD3D","CD3E"))

FilterPlots(sample_name, G_RNA_low = 500, G_RNA_hi= Inf, U_RNA_low = 0, U_RNA_hi = 37000, percent_mt_hi = 10, percent_ribo_hi= 100, percent_disso_hi = 100)

CalculateDoublets(sample_name, dbs_thr='none', dbs_remove=TRUE)

Normalize(sample_name, variable_genes=2000)

CCScore(sample_name)

ApplyRegression(UMI= sample_name, variables= c("nCount_RNA", "percent_mt", "S.Score", "G2M.Score"), explore_PC=TRUE)

CalculateCluster(sample_name, 12, "human")

MakeAnnotation(sample_name, organism="human")
}


\name{CalculateDoublets}
\alias{CalculateDoublets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CalculateDoublets}
\description{
Calculate and highlights doublets in the provided dataset. Generate also plots to show distribution and percentage of them.
}
\usage{
CalculateDoublets(UMI, dbs_thr='none', dbs_remove=TRUE, out_folder=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{Input UMI object generated via Seurat}
  \item{dbs_thr}{Doublets threshold. Default is auto calculated}
  \item{dbs_remove}{Remove doublets if TRUE. Default is TRUE}
  \item{out_folder}{Output folder. Default is the working directory}
}
\details{
CalculateDoublets by default returns two graphs in the "01.QC_Plots" folder:

the first (01h) allows the visualization of the cell doublets and of the doublet scores on a UMAP projection.

the second (01i) returns histograms of observed and simulated doublet score distributions and the value of the automatically detected threshold for the simulated doublet score.

}
\value{
Returns a Seurat Object after doublets calculation.

If "dbs_remove = TRUE", cells labeled as doublets are removed from the object.
}
\references{
The main function used to identify doublets has been modified from the R code provided at https://rdrr.io/github/ChengxiangQiu/rscrublet/ by Chengxiang Qiu \email{<cxqiu@uw.edu>} based on the method implemented in the python module Scrublet

Samuel L. Wolock, Romain Lopez, Allon M. Klein, "Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data" (2019). \emph{Cell System}, Volume \bold{8},--281-291--.e9,ISSN 2405-4712,\url{https://doi.org/10.1016/j.cels.2018.11.005.}
}
\author{Jimmy Caroli, Francesco Grandi, Chengxiang Qiu}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## first run:
umi_object <- CalculateDoublets(umi_object, dbs_thr='none', dbs_remove=FALSE)

## second run:
CalculateDoublets(umi_object, dbs_thr=0.22, dbs_remove=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

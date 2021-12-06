\name{FilterPlots}
\alias{FilterPlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{FilterPlots}
\description{
Apply filters to input data. Creates plot with provided thresholds as input.
}
\usage{
FilterPlots(UMI, G_RNA_low = 0, G_RNA_hi = Inf, U_RNA_low = 0, U_RNA_hi = Inf, percent_mt_hi = 100, percent_ribo_hi = 100, percent_ribo_hi = 100, out_folder=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{Input UMI object generated via Seurat}
  \item{G_RNA_low}{Minimum number of genes detected in each cell. Default is 0}
  \item{G_RNA_hi}{Maximum number of genes detected in each cell. Default is Infinite}
  \item{U_RNA_low}{Minimum number of molecules detected within a cell. Default is 0}
  \item{U_RNA_hi}{Maximum number of molecules detected within a cell. Default is Infinite}
  \item{percent_mt}{Maximum percentage of mitochondrial genes. Default is 100}
  \item{percent_ribo}{Maximum percentage of ribosomal genes. Default is 100}
  \item{percent_disso}{Maximum percentage of dissociation genes. Default is 100}
  \item{out_folder}{Output folder. Default is working directory}
}
\details{
FilterPlots by default returns two graphs in the "01.QC_Plots" folder:

the first plot (01f) displays data distributions after filtering; filtering thresholds are displayed as gray lines.

the second graph (01g) shows data correlations, that highlight with colors, cells passing the filtering procedure; filtering thresholds are displayed as gray lines.
}
\value{
Returns a Seurat Object filtered accordingly to the user-defined set of thresholds.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Jimmy Caroli, Francesco Grandi}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
FilterPlots(umi_object)

FilterPlots(umi_object, G_RNA_low = 500, G_RNA_hi = 10000, U_RNA_low = 1000, U_RNA_hi = 70000, percent_mt_hi = 20, percent_ribo_hi = 50, percent_disso_hi = 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

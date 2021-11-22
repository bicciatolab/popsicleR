\name{CalculateCluster}
\alias{CalculateCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CalculateCluster
%%  ~~function to do ... ~~
}
\description{
Calculate cluster of cells based on immune markers. Generate several plots.
}
\usage{
CalculateCluster(UMI, dim_pca, organism, marker.list='none', PCA=TRUE, cluster_res=0.8, out_folder=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{Input UMI object generated via Seurat}
  \item{dim_pca}{Number of principal component to take into account.}
  \item{organism}{Input organism. Defines the immune marker list. Default is human. Can be mouse}
  \item{marker.list}{List of markers. Default is a immune marker list defined by the selected organism. Can be fed personal list}
  \item{PCA}{If TRUE calculates also the PCA. Default is TRUE}
  \item{cluster_res}{Cluster resolution. Default is 0.8}
  \item{out_folder}{Output folder. Default is the working directory}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Jimmy Caroli}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
CalculateCluster(umi_object, 30, "human")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

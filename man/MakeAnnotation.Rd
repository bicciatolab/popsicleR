\name{MakeAnnotation}
\alias{MakeAnnotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MakeAnnotation
%%  ~~function to do ... ~~
}
\description{
Annotate the input umi object using Human Primary Cell Atlas and Blueprint Encode from SingleR (human).
Annotate the input umi object using scMCA and ImmGen from SingleR (mouse).
Function generates PCA, UMAP and TSNE plots
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MakeAnnotation(UMI, organism=c("human","mouse"), marker.list='none', thresh=20, cluster_res=NULL, out_folder=getwd())
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{ The input umi object generated via Seurat}
  \item{organism}{ The organism on which perfom the analysis. Can be human or mouse.}
  \item{marker.list}{ Markers list used for the plots. Must be the same used in CalculateCluster. Default uses proprietary lists.}
  \item{thresh}{ Rare population number threshold. Default is 20}
  \item{cluster_res}{ Defines on which cluster resolution perform the annotation. Must be a clustering resolution already calculated by CalculateCluster in the previous step. Default is NULL.}
  \item{out_folder}{ Output folder. Default is the working directory}
}
\details{
MakeAnnotation by default returns several graphs in the "04.Annotation" dedicated subfolder:

graphs 04a, show single-cell annotations for each reference on reduced spaces projections (t-SNE and UMAP). Cells are colored according to the assigned cell type.

graphs 04b, show cluster annotations for each reference on reduced spaces projections (t-SNE and UMAP). Cells are colored according to the assigned cell type.

graphs 04c, report single-cell annotations for each reference on reduced spaces projections (t-SNE and UMAP). Cells are colored in red according to the assigned cell type. Each cell type is displayed individually.

graphs 04d, return dot plots reporting, in the cell types assigned using the different references, the markers genelist expression.

graphs 04e, are pie-chart plots displaying, for each cluster, the proportion of cells annotated with a specific label based on the different references.
}
\value{
Returns a Seurat Object after annotation.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Jimmy Caroli. Francesco Grandi}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
MakeAnnotation(umi_object, organism="human", cluster_res=0.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{MakeAnnotation}
\alias{MakeAnnotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MakeAnnotation
%%  ~~function to do ... ~~
}
\description{
MakeAnnotation annotates cells using SingleR and scMCA built-in references.

For human samples, annotation is performed basing on the Human Primary Cell Atlas (HPCA) and the Blueprint Encode (BpEn) references from celldex and SingleR packages.

For mouse samples, annotation is performed according to mouse bulk RNA-seq datasets (Mouse RNA-seq) and Immunological Genome Project (ImmGen) references from SingleR, and the single-cell Mouse Cell Atlas (scMCA) reference.

Annotations are displayed on PCA, UMAP and TSNE embeddings, and evaluated through further graphs.
}
\usage{
MakeAnnotation(UMI, organism=c("human","mouse"), marker.list='none', thresh=20, cluster_res=NULL, out_folder=getwd())
}



%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{ The input umi object generated via CalculateCluster function}
  \item{organism}{ The organism on which perfom the analysis. Can be human or mouse.}
  \item{marker.list}{ Markers list used for the plots. Must be the same used in CalculateCluster. Default uses proprietary lists.}
  \item{thresh}{ Rare population number threshold. Default is 20}
  \item{cluster_res}{ Defines on which cluster resolution perform the annotation. Must be a clustering resolution already calculated by CalculateCluster in the previous step. Default is NULL.}
  \item{out_folder}{ Output folder. Default is the working directory}
}
\details{
MakeAnnotation returns several graphs in the "04.Annotation" dedicated subfolder:

graphs 04a, show single-cell annotations for each reference on reduced spaces projections (t-SNE and UMAP). Cells are colored according to the assigned cell type.

graphs 04b, show cluster annotations for SingleR references on reduced spaces projections (t-SNE and UMAP). Cells are colored according to the clustering and labeled basing on the assigned cell type.

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
MakeAnnotation(umi_object, organism="human", marker.list= "markers.txt", cluster_res=0.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

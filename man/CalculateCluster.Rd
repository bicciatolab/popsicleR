\name{CalculateCluster}
\alias{CalculateCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CalculateCluster
%%  ~~function to do ... ~~
}
\description{
Calculate cluster of cells based on immune markers. Clusters will be displayed across tSNE and UMAP plots.
}
\usage{
CalculateCluster(UMI, dim_pca, organism=c("human","mouse"), marker.list='none', PCA=TRUE, cluster_res=0.8, out_folder=getwd())
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{Input UMI object generated via Seurat}
  \item{dim_pca}{Number of principal component to take into account.}
  \item{organism}{Input organism to define the marker list to use. Can be human or mouse}
  \item{marker.list}{List of markers. Default is 'none' and a proprietary marker list is applied for the selected organism. User can fed a custom genelist specifying the file_name of a tab delimited file stored in the working directory.}
  \item{PCA}{If TRUE calculates the PCA. Default is TRUE}
  \item{cluster_res}{Cluster resolution. Default is 0.8. User can specify, through a vector, multiple clustering resolution to explore. When input more than one resolution, time consuming steps (e.g. Seurat::FindAllMarkers) are not performed.}
  \item{out_folder}{Output folder. Default is the working directory}
}
\details{
CalculateCluster by default returns several graphs in the "03.Clustering" dedicated subfolder:

graphs 03a, show for each user-provided clustering resolution projections on reduced spaces (PCA, t-SNE, and UMAP) colored according to the clustering.

graph 03b, computed only when more clustering resolution are provided, displays a clustering tree that shows how cells are assigned to the clusters at the various clustering resolutions

graph 03c, reports a phylogenetic tree showing correlation between the identified clusters for each clustering resolution.

graphs 03d, explore data distribution based on QC features, MALAT1 and GAPDH expression and doublets estimation through projections on reduced spaces (PCA, t-SNE, and UMAP) final embeddings.

graph 03e, computed only when a single clustering resolution is provided by the user, returns violin plots reporting the expression of the 2 top marker genes for each cluster.

graphs 03f, computed only when a single clustering resolution is provided by the user, report the expression of the 2 top marker genes for each cluster on reduced spaces projections (t-SNE and UMAP).

graph 03g, computed only when a single clustering resolution is provided by the user, plots a heatmap that resume for each cluster the top10-markers expression values.

graphs 03h, show marker genes expression on reduced spaces projections (t-SNE and UMAP).

graph 03i, returns violin plots of marker genes expression for each clustering resolution provided.

graph 03j, displays dot plots reporting marker genes expression in the identified clusters for each user-provided clustering resolution.
}
\value{
Returns a Seurat Object after clustering.
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
CalculateCluster(umi_object, 12, organism= "human", cluster_res=c(0.4, 0.6, 0.8))

CalculateCluster(umi_object, 12, organism= "human", cluster_res= 0.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

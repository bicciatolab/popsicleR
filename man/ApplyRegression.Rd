\name{ApplyRegression}
\alias{ApplyRegression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ApplyRegression
%%  ~~function to do ... ~~
}
\description{
Function to perform the regression. Input will require variable/s on which the regression will be calculated
}
\usage{
ApplyRegression(UMI, organism=c("human","mouse"), variables='none', explore_PC=FALSE, out_folder=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{Input UMI object, generated via Seurat}
  \item{organism}{The organism on which perfom the analysis. Can be human or mouse. Organism must be specified to calculate the Cell Cycle Score based on a tailored proprietary genelist}
  \item{variables}{List (or single digit) of variables to be regressed on. Variables can be: nFeature_RNA, nCount_RNA, percent_mt, S.Score, G2M.Score. Default is 'none' as to not apply any regression}
  \item{explore_PC}{When TRUE performs additional investigations and generates more graphs. Default is FALSE.}
  \item{out_folder}{Output folder. Default is the working directory}
}
\details{
ApplyRegression by default returns several graphs in the "02.PreProcessing" dedicated subfolder:

graph 02c, when variables on which perform regression are not specified returns, in the "No_Regression" folder, projections on reduced spaces (PCA, t-SNE, and UMAP using 20 principal components) to observe data distribution without regression.
When variables to regress are specified returns projections on reduced spaces to inspect the effect of variable regression.

graph 02d displays top genes associated with principal components.

graph 02e shows a heatmap focusing on a principal component. Both cells and genes are sorted by their principal component scores.

graph 02f returns the results of the JackStraw analysis for PCA significance.

graph 02g reports the standard deviations of the principle components for easy identification of an elbow in the graph. This elbow often corresponds well with the significant PCs and is much faster to run.
}
\value{
Returns a Seurat Object after regression.
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
ApplyRegression(UMI= umi_object, organism= "human", variables= 'none', explore_PC=FALSE)

ApplyRegression(UMI= umi_object, organism= "human", variables= c("nCount_RNA", "percent_mt", "S.Score", "G2M.Score"), explore_PC=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

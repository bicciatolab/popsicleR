\name{ApplyRegression}
\alias{ApplyRegression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ApplyRegression
%%  ~~function to do ... ~~
}
\description{
Perform data scaling and optionally regression on normalized data. User can specify variable/s on which calculate the regression.
Before scaling cell cycle scores are calculated.  When regression is required, this function should be run twice sequentially, first time to explore variables distribution and subsequently to regress on appropriate variables.
}
\usage{
ApplyRegression(UMI, organism=c("human","mouse"), variables='none', explore_PC=FALSE, out_folder=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{Input UMI object, generated via Normalize function}
  \item{organism}{The organism on which perfom the analysis. Can be human or mouse.}
  \item{variables}{String or character vector specifying variable/s on which perform regression. Commonly, variables used for regression are: nFeature_RNA, nCount_RNA, percent_mt, S.Score, G2M.Score. Default is 'none' as to not apply any regression}
  \item{explore_PC}{When TRUE performs additional investigations and generates further exploration graphs. Default is FALSE.}
  \item{out_folder}{Output folder. Default is the working directory}
}
\details{
ApplyRegression returns several graphs in dedicated subfolders inside the "02.PreProcessing" directory. When no regression variables are specified, a "No_Regression" subfolder is generated. Otherwise, a dedicated subfolder is generated for each combination of regression variables.

graph 02c, returns projections on reduced spaces (PCA, t-SNE, and UMAP using 20 principal components) to observe variables distribution.

Furthermore, when explore_PC=TRUE graphs from 02d to 02g are generated in order to choose the optimal number of principal components to consider for further analysis.

graph 02d displays top genes associated with the first two principal components.

graph 02e shows a heatmap focusing on a principal component. Both cells and genes are sorted by their loadings for each principal component.

graph 02f returns the results of the JackStraw analysis for PCA significance.

graph 02g reports the standard deviation observed for each principal component.
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

ApplyRegression(UMI= umi_object, organism= "human", variables= c("nFeature_RNA","nCount_RNA", "percent_mt", "S.Score", "G2M.Score"), explore_PC=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

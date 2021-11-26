\name{Annotate}
\alias{Annotate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Annotate
%%  ~~function to do ... ~~
}
\description{
Annotate the input umi object using Human Primary Cell Atlas and Blueprint Encode (from SingleR) if the organism selected is human.
Annotate the input umi object using scMCA and ImmGen (from SingleR) if the organism selected is mouse.
Function generate PCA, UMAP and TSNE plots
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Annotate(umi_scaled, organism, thresh=20, out_folder=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{umi_scaled}{ The input umi object generated via Seurat}
  \item{organism}{ Organism for the annotation. Select human or mouse}
  \item{thresh}{ Rare population number threshold. Default is 20}
  \item{out_folder}{ Output folder. Default is the working directory}
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
Annotate(UMI, 'human', 30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

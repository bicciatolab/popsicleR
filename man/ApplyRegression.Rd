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
ApplyRegression((UMI, variables='none', explore_PC=FALSE,  out_folder=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{Input UMI object, generated via Seurat}
  \item{variables}{List (or single digit) of variables to be regressed on. Variables can be: nFeature_RNA, nCount_RNA, percent_mt, S.Score, G2M.Score. Default is 'none' as to not apply any regression}
  \item{explore_PC}{When TRUE performs additional investigations and generates more graphs. Default is FALSE.}
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
ApplyRegression(umi_object, 'nFeature_RNA')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

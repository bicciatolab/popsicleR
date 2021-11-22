\name{CalculateDoublets}
\alias{CalculateDoublets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CalculateDoublets}
\description{
Calculate doublets using external python script
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
CalculateDoublets(umi_object)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

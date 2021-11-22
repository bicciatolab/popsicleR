\name{CCScore}
\alias{CCScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CCScore}
\description{
Calculate the Cell Cycle Score
}
\usage{
CCScore(UMI, organism='human')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{Input UMI object generated via Seurat}
  \item{organism}{Organism. Default is human. Mouse is available}
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
CCScore(umi_object)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

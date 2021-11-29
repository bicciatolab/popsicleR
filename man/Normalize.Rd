\name{Normalize}
\alias{Normalize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Normalize}
\description{
Perform data normalization of the provided input. User must define the number of variable genes to be investigated
after normalization.
}
\usage{
Normalize(UMI, variable_genes=2000, out_folder=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{UMI}{Input UMI object generated via Seurat}
  \item{variable_genes}{Number of highly variable genes investigated after normalization. Default is 2000}
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
\author{Jimmy Caroli, Francesco Grandi}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
seurat_object <- Normalize(seurat_object, variable_genes=2000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

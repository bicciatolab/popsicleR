\name{plotGene}
\alias{plotGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plotGene
%%  ~~function to do ... ~~
}
\description{
Creates feature plots for the selected input gene, starting from the input umi object.
}
\usage{
plotGene(genelist, umi, dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genelist}{Genelist. Will be generated plots for genes in genelist.}
  \item{umi}{Input UMI object obtained from Seurat.}
  \item{dir}{Output directory for generated plots.}
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
plotGene(c('PTPRC','TP53'), UMI, "/output_directory")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

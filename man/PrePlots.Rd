\name{PrePlots}
\alias{PrePlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PrePlots}
\description{
Perfoms preprossing analyses and quality controls on the provided dataset.
Reports information on the dataset and outputs plots describing quality of the data to be further analyzed. 
}
\usage{
PrePlots(sample_name, input_data, genelist, percentage=0.1, gene_filter=200, cellranger=TRUE, input_matrix=NULL, organism="human", out_folder=getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_name}{Name of the sample investigated}
  \item{input_data}{Input folder location}
  \item{genelist}{List of gene for singular gene plots}
  \item{percentage}{Percentage of cell kept in the analysis. Default is 0.1}
  \item{gene_filter}{Number of genes to keep in the analysis. Default is 200}
  \item{cellranger}{Data were generated via cellranger. Default is TRUE}
  \item{input_matrix}{Input matrix data. Not generated via cellranger. Default is NULL}
  \item{organism}{Default is human. Can be mouse.}
  \item{out_folder}{Output folder. Default is working directory}
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
PrePlots('breast_single_cell', TRUE, path_to_data, 0.1, 200, c('TP53','PTEN'))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

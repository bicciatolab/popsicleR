\name{PrePlots}
\alias{PrePlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PrePlots}
\description{
Perfoms preprossing analyses and quality controls on the provided dataset.
Reports information on the dataset and outputs plots describing quality of the data to be further analyzed.
}
\usage{
PrePlots(sample_name, input_data, genelist=NULL, percentage=0.1, gene_filter=200, cellranger=TRUE, organism=c("human","mouse"), out_folder=getwd())
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_name}{Name of the sample investigated}
  \item{input_data}{Input data, if cellranger=TRUE must be a path to a folder containing cellranger output. Else, if input data are not generated via cellranger, it must be a .txt matrix with gene names in the first column and cell names in the first row.}
  \item{genelist}{List of gene for singular gene plots. By default is NULL and a proprietary genelist is used to generate plots. A personal genelist can be fed as characters vector}
  \item{percentage}{Percentage of cell kept in the analysis. Default is 0.1}
  \item{gene_filter}{Number of genes to keep in the analysis. Default is 200}
  \item{cellranger}{Specify if data were generated via cellranger. Default is TRUE}
  \item{organism}{The organism on which perfom the analysis. Can be human or mouse.}
  \item{out_folder}{Output folder. Default is working directory}
}
\details{
PrePlots by default returns several graphs in the "01.QC_Plots" folder:

violin (1a), density (1b), and scatter (1c) plots on number of genes (nFeature_RNA), number of UMI (nCount_RNA), the percentage of reads that map to the mitochondrial genome (percent_mt), to ribosomal (percent_ribo) and dissociation genes (percent_disso).

PrePlots also returns the distributions (1d) and the cell abundance (1e) of specific cell populations identified by a list of user-defined marker genes.
}
\value{
Returns a Seurat Object generated from input data.

Genes expressed in less than 10 percent of cells and cells that express less than 200 genes are filtered-out from input data.
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
PrePlots('breast_single_cell', input_data="/path/to/cellranger_output_directory/", genelist=NULL, 0.1, 200, cellranger = TRUE, c('TP53','PTEN'), organism="human")

PrePlots("mouse_sample", input_data="/path/to/matrix.txt", cellranger=FALSE, organism= "mouse")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

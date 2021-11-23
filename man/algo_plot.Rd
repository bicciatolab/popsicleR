\name{algo_plot}
\alias{algo_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{algo_plot
%%  ~~function to do ... ~~
}
\description{
Creates plots using the algorithm defined as input.
}
\usage{
algo_plot(dir, data, type, plot_name, res)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{Directory in which the plots will be saved.}
  \item{data}{Input UMI object.}
  \item{type}{Reduction type.}
  \item{plot_name}{Name of the plot}
  \item{res}{Resolution}
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
\author{Jimmy Caroli
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
algo_plot(Cluster_dir, umi_scaled, "umap", "04a_clusters", 0.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{quantbox}
\alias{quantbox}
\title{
  Redefine Boxplot to Show Quantiles
}
\description{
  Hijack boxplot code and replace \code{boxplot.stats} with \code{myboxplot.stats}.
}
\usage{
quantbox(x, use.cols = TRUE, ...)
}
\arguments{
  \item{x}{\code{numeric} -- data from which quantile plots are created;
    either a numeric vector, a numeric matrix, or a single list containing numeric vectors}
  \item{use.cols}{\code{logical} -- if \code{TRUE}, plot columns, otherwise, plot rows}
  \item{\dots}{\code{dots} -- further arguments to \code{\link[graphics]{boxplot}}.}
}
\details{
  The objects \code{'myboxplot.stats'} and \code{'qboxplot'} (boxplot modified to render quantiles)
  are located in the package environment: \code{.PBSmodEnv}. The function \code{'quantbox'}
  grabs \code{'qboxplot'} using the accessor function \code{'tget'}.
}
\value{
  Invisibly returns a list object containing:
  \tabular{ll}{
    \code{stats} \tab matrix of quantiles,\cr
    \code{n}     \tab a vector with the number of observations in each group,\cr
    \code{conf}  \tab a matrix where each column contains the lower and upper extremes of the notch;\cr
                 \tab \code{stats[3] + c(-1.58, 1.58) * diff(stats[c(2, 4)])/sqrt(n)}\cr
    \code{out}   \tab the values of any data points which lie beyond the extremes of the whiskers,\cr
    \code{group} \tab a vector of the same length as \code{out} whose elements indicate to which group the outlier belongs,\cr
    \code{names} \tab names for each box plotted.\cr
  }
}
\references{
  \href{https://stat.ethz.ch/pipermail/r-help/2011-April/275384.html}{R Help originally on Nabble [now defunct]}
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\note{
  This function is a clone of the \pkg{PBStools} function \code{'quantBox'}.
  Why would anyone choose a boxplot over a quantile plot?
}
\seealso{
In package \pkg{graphics}:\cr
  \code{\link[graphics]{boxplot}},
  \code{\link[graphics]{boxplot.matrix}}\cr
In package \pkg{PBSmodelling}:\cr
  \code{\link[PBSmodelling]{tget}}\cr
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{compB0}},
  \code{\link[PBSawatea]{panelBoxes}},
  \code{\link[PBSawatea]{plt.ageResidsPOP}},
  \code{\link[PBSawatea]{plt.cohortResids}},
  \code{\link[PBSawatea]{plt.yearResidsPOP}}
}
\examples{
\dontrun{
local(envir=.PBSmodEnv,expr={
  rien = list('1997'=rnorm(100,10),'2005'=rnorm(100,12),'2016'=rnorm(100,14))
  quantbox(rien,outline=FALSE)
})
}
}
\keyword{hplot}
\concept{plotFuns.r}


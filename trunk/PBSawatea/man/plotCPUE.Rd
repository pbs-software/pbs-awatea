\name{plotCPUE}
\alias{plotCPUE}
\title{
  Plot CPUE and Add Error Bars
}
\description{
  Plot CPUE and fit with error bars.
}
\usage{
plotCPUE(obj, main="", save=NULL, bar=1.96, yLim=NULL, ...)
}
\arguments{
  \item{obj}{data frame of CPUE indices from Awatea's results file
    (\emph{e.g.}, \code{currentRes$CPUE}).}
  \item{main}{title for figure}
  \item{save}{\bold{not used}}
  \item{bar}{standard deviation of the normal distribution 
    (1.96 is the approximate value of the 97.5 percentile point).}
  \item{yLim}{limits of the y-axis.}
  \item{\dots}{\bold{not used}}
}
\value{
  A postscript file: \cr
  \code{CPUEser.eps...}{CPUE indices with error bars.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Copied code from \code{plotIndexNotLattice}.
}
\seealso{
  \code{\link{plotIndexNotLattice}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

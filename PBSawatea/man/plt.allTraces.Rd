\name{plt.allTraces}
\alias{plt.allTraces}
\title{
  Plot MCMC Traces
}
\description{
  Plot traces from MCMC samples.
}
\usage{
plt.allTraces(obj, bioYrList=NULL, recYrList=NULL, save=TRUE)
}
\arguments{
  \item{obj}{\code{numeric} -- vector of MCMC samples.}
  \item{bioYrList}{\code{numeric} -- years to plot spawning biomass traces.}
  \item{recYrList}{\code{numeric} -- years to plot recruitment traces.}
  \item{save}{\code{logical} -- if \code{TRUE}, save figure to a raster file (\code{.jpg}).}
}
\author{
First created by:\cr
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew M. Edwards}, PBS, DFO, Nanaimo BC

Maintained by:\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\note{
  Appears to be some figure requested by PJS.
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{panelChains}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

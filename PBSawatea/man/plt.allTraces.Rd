\name{plt.allTraces}
\alias{plt.allTraces}
\title{
  Plot MCMC Traces
}
\description{
  Plot  traces from MCMC samples.
}
\usage{
plt.allTraces(obj, bioYrList=NULL, recYrList=NULL, save=TRUE)
}
\arguments{
  \item{obj}{vector of MCMC samples.}
  \item{bioYrList}{years to plot spawning biomass traces.}
  \item{recYrList}{years to plot recruitment traces.}
  \item{save}{logical: if \code{TRUE}, save figure to a raster file (\code{.jpg}).}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Appears to be some figure requested by PJS.
}
\seealso{
  \code{\link{plotChains}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

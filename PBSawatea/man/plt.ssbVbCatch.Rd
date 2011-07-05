\name{plt.ssbVbCatch}
\alias{plt.ssbVbCatch}
\title{
  Plot Annual Spawning and Vulnerable Biomass
}
\description{
  Plot MPD values of spawning biomass (SB) and vulnerable biomass(VB), as well as catch, against year.
}
\usage{
plt.ssbVbCatch(obj, x1=1966, xLim=c(1954,2005), yLim=c(0,25000))
}
\arguments{
  \item{obj}{an Awatea results object from \code{load.allResFiles}.}
  \item{x1}{year to start plotting SB and VB lines.}
  \item{xLim}{limits of the x-axis.}
  \item{yLim}{limits of the y-axis.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  This analysis uses the MPD (mode of the posterior distribution) values for \eqn{B} and \eqn{V}.
}
\seealso{
  \code{\link{load.allResFiles}}, \code{\link{get.resFile}}, \code{\link{plt.expRate}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

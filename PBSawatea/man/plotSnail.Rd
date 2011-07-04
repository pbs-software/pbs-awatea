\name{plotSnail}
\alias{plotSnail}
\title{
  Plot Snail Trails of Exploitation vs. Biomass
}
\description{
  Plot the historical progression of the ratio \eqn{u_t/u_{MSY}}{ut/umsy}
  against \eqn{B_t/B_{MSY}}{Bt/Bmsy}.
}
\usage{
plotSnail(BoverBmsy, UoverUmsy, p=c(0.1,0.9), 
          xLim=NULL, yLim=NULL, Lwd=2)
}
\arguments{
  \item{BoverBmsy}{numeric matrix of \eqn{B_t}{Bt} over \eqn{B_{MSY}}{Bmsy}.}
  \item{UoverUmsy}{numeric matrix of \eqn{u_t}{ut} over \eqn{u_{MSY}}{umsy}.}
  \item{p}{quantiles to show the bulk of the distribution.}
  \item{xLim}{limits of the x-axis.}
  \item{yLim}{limits of the y-axis.}
  \item{Lwd}{line width of the snail trail.}
}
\details{
  The graph attempts to show the time history of the exploitation rate 
  compared to the spawning biomass using a precautionary framework recast
  in Bayesian terms.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  The term \emph{snail trail} comes from PJS.
}
\seealso{
  \code{\link{plotBVBnorm}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

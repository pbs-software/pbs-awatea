\name{refPoints}
\alias{refPoints}
\alias{refPointsB0}
\title{
  Calculate Reference Points
}
\description{
  Calculate reference points relative to either \eqn{B_{MSY}}{Bmsy} or \eqn{B_0}{B0}.
}
\usage{
refPoints(mcmcObj=currentMCMC, projObj=currentProj, 
          msyObj=currentMSY, refLevels=c(0.4,0.8,1))

refPointsB0(mcmcObj=currentMCMC, projObj=currentProj,
          B0Obj=B0.MCMC, refLevels=B0refLevels, refNames=B0refNames)
}
\arguments{
  \item{mcmcObj}{MCMC list object (\emph{e.g.}, \code{currentMCMC}).}
  \item{projObj}{projected biomass list object (\emph{e.g.}, \code{currentProj}).}
  \item{msyObj}{MSY list object (\emph{e.g.}, \code{currentMSY}).}
  \item{refLevels}{reference levels relative to \eqn{B_{MSY}}{Bmsy} (or \eqn{B_0}{B0}).}
  \item{B0Obj}{vector of \eqn{B_0}{B0} MCMC values (\emph{e.g.}, \code{B0.MCMC}).}
  \item{refNames}{names of the \eqn{B_0}{B0} reference levels \code{refLevels}.}
}
\value{
  List of reference points relative to either \eqn{B_{MSY}}{Bmsy} or \eqn{B_0}{B0}.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Call from \code{Sweave} as  \code{refPoints()} or, in full: \cr
  \code{refPoints(currentMCMC, currentProj, currentMSY, refLevels=c(0.4,0.8,1))}
}
\seealso{
  \code{\link{calc.refVal}}
}
\keyword{arith}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

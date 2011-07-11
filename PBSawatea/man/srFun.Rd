\name{srFun}
\alias{srFun}
\title{
  Stock Recruitment Function
}
\description{
  Take a vector of spawners in year \eqn{t}-1 and calculate 
  recruits in year \eqn{t}. 
}
\usage{
srFun(spawners, h = h.mpd, R0 = R0.mpd, B0 = B0.mpd)
}
\arguments{
  \item{spawners}{a vector of spawners where either : \cr
     each element corresponds to spawners in year \eqn{t}-1 or \cr
     the vector calculates a single year but multiple MCMCs.}
  \item{h}{steepness parameter value.}
  \item{R0}{recruitment at \eqn{t} = 0 (virgin conditions).}
  \item{B0}{spawning biomasss at virgin conditons.}
}
\details{
  (AME wording) To input a vector of spawners in year \eqn{t}-1
  and calculate recruits in year \eqn{t}. Output for recruits
  is vector, each element corresponds to spawners the
  the year before, so will usually want to shift the
  output by 1 so that recruits in year t are based on
  spawners in year \eqn{t}-1.

  Can also have each input as a vector (used when
  calculating a single year but multiple MCMCs, as in
  first year of projections is based on penultimate year
  of MCMC calcualtions).
}
\value{
  A vector of recruitments in year \eqn{t}.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  This function was originally a subfunction in \code{plt.mpdGraphs}.
}
\seealso{
  \code{\link{plt.mpdGraphs}}
}
\keyword{ts}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


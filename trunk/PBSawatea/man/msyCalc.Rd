\name{msyCalc}
\alias{msyCalc}
\title{
  Calculate the Maximum Sustainable Yield
}
\description{
  Load in \code{MSY.out} and calculate the MSY (maximum sustainable yield).
}
\usage{
msyCalc(dir = getwd(), error.rep = 1)
}
\arguments{
  \item{dir}{working directory.}
  \item{error.rep}{numeric: if \code{1}, report errors (reaching bounds), if \code{0} do not.}
}
\value{
  Returns a list object containing: \cr
    \code{yield...}{maximum sustainable yield,} \cr
    \code{u.......}{exploitation rate at MSY,} \cr
    \code{VB......}{vulnerable biomass at MSY,} \cr 
    \code{B.......}{spawning biomass at MSY,} \cr
    \code{nProj...}{numnber of projections needed to reach MSY.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  See \code{msyTestCreating.r} for full details when figuring this out.
}
\seealso{
  \code{\link{findTarget}}
}
\keyword{interation}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

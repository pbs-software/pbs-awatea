\name{compB0}
\alias{compB0}
\title{
  Compare Reference Criteria / Points in B0 Space
}
\description{
  Compare COSWEIC reference criteria and DFO reference points in terms of \eqn{B_0}{B0}. \cr
  The figure concept comes from Chris Woods (PBS).
}
\usage{
compB0(B, Mnam=c("Est M", "Fix M"), ratios=c(0.4, 0.8), 
    incl.Bt=TRUE, boxwidth=0.6, figgy=FALSE, width=12, height=9)
}
\arguments{
  \item{B}{list of list of MCMC samples (see \bold{Details}); 
    the first level of the list is the model run, while the second level contains 
    MCMC samples (one of which should be \eqn{B_0}{B0} which acts as the divisor to the other MCMCs).}
  \item{Mnam}{model names for the boxplot.}
  \item{ratios}{reference levels of \eqn{B_{MSY}}{Bmsy} (usually 0.4 and 0.8).}
  \item{incl.Bt}{logical: if \code{TRUE}, include the current spawning stock biomass in terms of \eqn{B_0}{B0}.}
  \item{boxwidth}{width of the box in x-units.}
  \item{figgy}{logical: if \code{TRUE}, send figure to three output files (\code{.eps},\code{.png}, and \code{.wmf}).}
  \item{width}{width of the output files in inches.}
  \item{height}{height of the output files in inches.}
}
\details{
  An example of the input list \code{B}: \cr
  \code{List of 2} \cr
  \code{..$ 29.01:List of 3} \cr
  \code{....$ B0.MCMC  : num [1:1000]} \cr
  \code{....$ Bt.MCMC  : num [1:1000]} \cr
  \code{....$ Bmsy.MCMC: num [1:1000]} \cr
  \code{..$ 30.01:List of 3} \cr
  \code{....$ B0.MCMC  : num [1:1000]} \cr
  \code{....$ Bt.MCMC  : num [1:1000]} \cr
  \code{....$ Bmsy.MCMC: num [1:1000]}
}
\value{
  Invisibly returns the list object \code{xBox} used to create the boxplot.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Uses a modified version of \code{boxplot} called \code{\link[PBSawatea]{plotBox}}.
}
\seealso{
  \code{\link[PBSawatea]{plotBox}}, \code{\link[PBSawatea]{importMCMC}}, \code{\link[PBSawatea]{msyCalc}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

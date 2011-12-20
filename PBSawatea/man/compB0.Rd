\name{compB0}
\alias{compB0}
\title{
  Compare Reference Criteria and Points Relative to B0
}
\description{
  Compare COSWEIC reference criteria and DFO reference points relative to \eqn{B_0}{B0}. \cr
  The figure concept comes from Chris Woods (PBS).
}
\usage{
compB0(B, Mnams=NULL, ratios=c(0.4,0.8), 
   include=list(A1=TRUE, A2=TRUE, SSPM=TRUE, Bmsy=TRUE, Bt=TRUE),
   t.yr=2011, boxwidth=0.6, figgy=FALSE, width=12, height=9, ...)
}
\arguments{
  \item{B}{list of list of MCMC samples (see \bold{Details}); 
    the first level of the list is the model run, while the second level contains 
    MCMC samples (one of which should be \eqn{B_0}{B0} which acts as the divisor to the other MCMCs).}
  \item{Mnams}{optional model names for the boxplot.}
  \item{ratios}{reference levels of \eqn{B_{MSY}}{Bmsy} (usually 0.4 and 0.8).}
  \item{include}{list of logicals specifying whether to plot 
    COSEWIC's reference criteria \code{A1} and \code{A2},
    the Schaefer surplus production model \code{SSPM},
    the biomass at maximum sustainable yield \code{Bmsy}, and
    the biomass of a given time period \code{Bt}, where \code{t.yr} specifies the year.
    All values are cast in terms of \eqn{B_0}{B0}.}
  \item{t.yr}{numeric year that represents the time period for \eqn{B_t}{Bt}.}
  \item{boxwidth}{width of the box in x-units.}
  \item{figgy}{logical: if \code{TRUE}, send figure to four output files (\code{.eps}, \code{.pdf}, \code{.png}, and \code{.wmf}).}
  \item{width}{width of the output files in inches.}
  \item{height}{height of the output files in inches.}
  \item{\dots}{additional values for \code{\link[PBSawatea]{plotBox}}'s argument \code{pars}.}
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

  The function creates a figure comparing COSEWIC criteria and reference points for 
  the Schaefer surplus production model and specified ratios of \eqn{B_{MSY}}{Bmsy}
  from catch-at-age models. The coordinate space (\emph{y}-axis) is relative to \eqn{B_0}{B0}.
}
\value{
  Invisibly returns a list object of \code{xBox} and \code{BarBox} used to create the boxplot.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Uses a modified version of \code{boxplot} called \code{\link[PBSawatea]{plotBox}}.
}
\seealso{
  \code{\link[PBSawatea]{plotBox}}, \code{\link[scapeMCMC]{importMCMC}}, \code{\link[PBSawatea]{msyCalc}} \cr
  \code{\link[grDevices]{plotmath}}, \code{\link[graphics]{boxplot}}, \code{\link[graphics]{bxp}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

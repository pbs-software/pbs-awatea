\name{compBmsy}
\alias{compBmsy}
\title{
  Compare Biomass Posteriors Relative to Bmsy
}
\description{
  Compare posterior distributions of biomass from MCMCs for 
  species and model runs.
}
\usage{
compBmsy(Bspp, spp="POP", Mnams=c("Est M","Fix M"), 
    ratios=c(0.4,0.8), t.yr=2011, figgy=FALSE, width=12, height=9, ...)
}
\arguments{
  \item{Bspp}{list of species MCMC results, which are lists of model runs, 
    each listing data frames of \eqn{B_t}{Bt} and \eqn{B_{MSY}}{Bmsy}.}
  \item{spp}{species code(s) (e.g., 3-letter), which are the first-level names of \code{Bspp}.}
  \item{Mnams}{optional model names for the boxplots.}
  \item{ratios}{reference levels of \eqn{B_{MSY}}{Bmsy} (usually 0.4 and 0.8).}
  \item{t.yr}{numeric year that represents the time period for \eqn{B_t}{Bt}.}
  \item{figgy}{logical: if \code{TRUE}, send figure to four output files (\code{.eps}, \code{.pdf}, \code{.png}, and \code{.wmf}).}
  \item{width}{width of the output files in inches.}
  \item{height}{height of the output files in inches.}
  \item{\dots}{additional values for \code{\link[PBSawatea]{plotBox}}'s argument \code{pars}.}
}
\details{
  An example of the input list \code{Bspp}: \cr
  \code{List of 2} \cr
  \code{..$ POP:List of 2} \cr
  \code{....$ run23:List of 3} \cr
  \code{......$ B0.MCMC  : num [1:1000]} \cr
  \code{......$ Bt.MCMC  : num [1:1000]} \cr
  \code{......$ Bmsy.MCMC: num [1:1000]} \cr
  \code{....$ run16:List of 3} \cr
  \code{......$ B0.MCMC  : num [1:1000]} \cr
  \code{......$ Bt.MCMC  : num [1:1000]} \cr
  \code{......$ Bmsy.MCMC: num [1:1000]} \cr
  \code{....- attr(*, "spp")= chr "POP"} \cr
  \code{..$ YMR:List of 2} \cr
  \code{....$ 29.01:List of 3} \cr
  \code{......$ B0.MCMC  : num [1:1000]} \cr
  \code{......$ Bt.MCMC  : num [1:1000]} \cr
  \code{......$ Bmsy.MCMC: num [1:1000]} \cr
  \code{....$ 30.01:List of 3} \cr
  \code{......$ B0.MCMC  : num [1:1000]} \cr
  \code{......$ Bt.MCMC  : num [1:1000]} \cr
  \code{......$ Bmsy.MCMC: num [1:1000]} \cr
  \code{....- attr(*, "spp")= chr "YMR"} \cr

  The function creates a set of horizontal boxes delimited by the 
  quantiles (0.025, 0.25, 0.5, 0.75, 0.975) that illustrate the posterior
  MCMC samples of biomass relative to \eqn{B_{MSY}}{Bmsy}. The default 
  reference points (0.4\eqn{B_{MSY}}{Bmsy} and 0.8\eqn{B_{MSY}}{Bmsy}) are 
  shown by vertical dashed lines.
}
\value{
  Invisibly returns the boxplot list object \code{Bmsy}.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC.
}
\seealso{
  \code{\link[PBSawatea]{compB0}}, \code{\link[PBSawatea]{plotBox}}
}
\keyword{hplot}


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
   t.yr=2011, boxwidth=0.6, figgy=FALSE, width=12, height=9,
   pngres=400, lang=c("e","f"), ...)
}
\arguments{
  \item{B}{\code{list} -- list of list of MCMC samples (see \bold{Details}); 
    the first level of the list is the model run, while the second level contains 
    MCMC samples (one of which should be \eqn{B_0}{B0} which acts as the divisor to the other MCMCs).}
  \item{Mnams}{\code{character} -- optional model names for the boxplot.}
  \item{ratios}{\code{numeric} -- reference levels of \eqn{B_{MSY}}{Bmsy} (usually 0.4 and 0.8).}
  \item{include}{\code{logical} -- list of logicals specifying whether to plot 
    COSEWIC's reference criteria \code{A1} and \code{A2},
    the Schaefer surplus production model \code{SSPM},
    the biomass at maximum sustainable yield \code{Bmsy}, and
    the biomass of a given time period \code{Bt}, where \code{t.yr} specifies the year.
    All values are cast in terms of \eqn{B_0}{B0}.}
  \item{t.yr}{\code{numeric} -- numeric year that represents the time period for \eqn{B_t}{Bt}.}
  \item{boxwidth}{\code{numeric} -- width of the box in x-units.}
  \item{figgy}{\code{logical} -- if \code{TRUE}, send figure to four output files (\code{.eps}, \code{.pdf}, \code{.png}, and \code{.wmf}).}
  \item{width}{\code{numeric} -- width of the output files in inches.}
  \item{height}{\code{numeric} -- height of the output files in inches.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) of the \code{png} file, if it's created.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
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
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-25}
}
\note{
  Uses a modified version of \code{boxplot} called \code{\link[PBSawatea]{plotBox}}.
}
\seealso{
  \code{\link[PBSawatea]{plotBox}}, \code{\link[scape]{importMCMC}}, \code{\link[PBSawatea]{msyCalc}} \cr
  \code{\link[grDevices]{plotmath}}, \code{\link[graphics]{boxplot}}, \code{\link[graphics]{bxp}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

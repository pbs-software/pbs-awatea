\name{plt.quantBio}
\alias{plt.quantBio}
\alias{plt.quantBioBB0}
\title{
  Plot Quantile Boxes of MCMC and Projected Biomass
}
\description{
  Plots MCMC and projected biomass as quantile boxes, the former in black, the latter in red.
}
\usage{
plt.quantBio(obj, projObj=NULL, policy=NULL, 
     p=tcall(quants5), xyType="lines", 
     lineType=c(3,2,1,2,3), refLines=NULL, 
     xLim=NULL, yLim=NULL, userPrompt=FALSE, 
     save=TRUE, yaxis.lab="Spawning biomass", lang="e")

plt.quantBioBB0(obj, projObj=NULL, policy=NULL, 
     p=tcall(quants5), xyType="lines", 
     lineType=c(3,2,1,2,3), refLines=NULL, 
     xLim=NULL, yLim=NULL, userPrompt=FALSE, save=TRUE, 
     main="", cex.main="", tcl.val=-0.2, 
     xaxis.by=1, yaxis.by=10000, xaxis.lab="Year", 
     yaxis.lab="Spawning biomass", lang="e")
}
\arguments{
  \item{obj}{an Awatea MCMC object (\emph{e.g.}, \code{currentMCMC}).}
  \item{projObj}{an Awatea projected biomass object (\emph{e.g.}, \code{currentProj}).}
  \item{policy}{numeric vector specifying catch policy.}
  \item{p}{quantiles to use from the biomass samples.}
  \item{xyType}{string specifying type of plot.}
  \item{lineType}{line types for the quantiles if \code{xyType="lines"}.}
  \item{refLines}{reference points.}
  \item{xLim}{limits of the x-axis.}
  \item{yLim}{limits of the y-axis.}
  \item{userPrompt}{logical: if \code{TRUE} prompts user before figure is drawn.}
  \item{save}{logical: if \code{TRUE} save figure as a raster file \code{.png}.}
  \item{main}{character string specifying a title for the plot.}
  \item{cex.main}{font size for figure title.}
  \item{tcl.val}{tick length.}
  \item{xaxis.by}{tick mark intervals for x-axis.}
  \item{yaxis.by}{tick mark intervals for y-axis.}
  \item{xaxis.lab}{label for x-axis.}
  \item{yaxis.lab}{label for y-axis.}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
}
\value{
  List of the reconstructed (MCMC) and projected results.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-26}
}
\note{
  \code{plt.quantBioBB0} performs similarly as for \code{plt.quantBio} but uses
  \eqn{B_t/B_0}{Bt/B0} instead of \eqn{B_t}{Bt}.
}
\seealso{
  \code{\link[PBSawatea]{plotVBcatch}},
  \code{\link[PBSawatea]{plotBVBnorm}}, \code{\link[PBSawatea]{plotRmcmcPOP}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

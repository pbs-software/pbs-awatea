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
   p=tcall(quants5), xyType="quantBox", lineType=c(3,2,1,2,3),
   delta=0.25, lwd=0.75, refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, main="", cex.main="",
   tcl.val=-0.1, xaxis.by=1, yaxis.by=10000,
   xaxis.lab="Year", yaxis.lab= "Spawning biomass", lang="e")
}
\arguments{
  \item{obj}{\code{list} -- an Awatea MCMC object (\emph{e.g.}, \code{currentMCMC}).}
  \item{projObj}{\code{list} -- an Awatea projected biomass object (\emph{e.g.}, \code{currentProj}).}
  \item{policy}{\code{numeric} -- vector specifying catch policy.}
  \item{p}{\code{numeric} -- quantiles to use from the biomass samples.}
  \item{xyType}{\code{character} -- string specifying type of plot.}
  \item{lineType}{\code{numeric} -- vector of line types for the quantiles if \code{xyType="lines"}.}
  \item{delta}{\code{numeric} -- x-offset for plotting boxes.}
  \item{lwd}{\code{numeric} -- line width device (see \code{'par'}).}
  \item{refLines}{\code{numeric} -- reference points.}
  \item{xLim}{\code{numeric} -- limits of the x-axis.}
  \item{yLim}{\code{numeric} -- limits of the y-axis.}
  \item{userPrompt}{\code{logical} -- if \code{TRUE} prompts user before figure is drawn.}
  \item{save}{\code{logical} -- if \code{TRUE} save figure as a raster file \code{.png}.}
  \item{main}{\code{character} -- string specifying a title for the plot.}
  \item{cex.main}{\code{numeric} -- font size for figure title.}
  \item{tcl.val}{\code{numeric} -- tick length.}
  \item{xaxis.by}{\code{numeric} -- tick mark intervals for x-axis.}
  \item{yaxis.by}{\code{numeric} -- tick mark intervals for y-axis.}
  \item{xaxis.lab}{\code{character} -- label for x-axis.}
  \item{yaxis.lab}{\code{character} -- label for y-axis.}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
}
\value{
  List of the reconstructed (MCMC) and projected results.
}
\author{
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew Edwards}, Research Scientist, PBS, DFO, Nanaimo BC\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2020-05-20}
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

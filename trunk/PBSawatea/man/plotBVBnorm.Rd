\name{plotBVBnorm}
%\alias{plotVBnorm}
\alias{plotBVBnorm}
\title{
  Plot Spawning and Vulnerable Biomass Relative to Virgin
}
\description{
  Plot spawning and vulnerable biomass boxplots relative to virgin levels 
  \eqn{B_0}{B0} and \eqn{V_0}{V0}, respectively.
}
\usage{
plotBVBnorm(mcmcObj, 
            p=c(0.025,0.25,0.5,0.75,0.975), 
            xyType="quantBox", 
            lineType=c(3,2,1,2,3), 
            refLines=NULL, xLim=NULL, yLim=NULL, 
            userPrompt=FALSE, save=TRUE, xLeg=0.7, yLeg=0.9, 
            yaxis.by=0.02, tcl.val=-0.2, 
            B.col="black", VB.col="black", ...)
}
%plotVBnorm( mcmcObj,
%            p=c(0.025,0.25,0.5,0.75,0.975),
%            xyType="quantBox",
%            lineType=c(3,2,1,2,3),
%            refLines=NULL, xLim=NULL, yLim=NULL,
%            userPrompt=FALSE, save=TRUE, xLeg=0.7, yLeg=0.9,
%            yaxis.by=0.02, tcl.val=-0.2,
%            B.col="black", VB.col="black", ...)

\arguments{
  \item{mcmcObj}{MCMC list object (\code{currentMCMC}).}
  \item{p}{quantiles to use in \code{quantBox}.}
  \item{xyType}{type of plot (currently only uses \code{quantBox}).}
  \item{lineType}{line types to use in \code{quantBox}.}
  \item{refLines}{reference lines to add to plot.}
  \item{xLim}{limits of the x-axis.}
  \item{yLim}{limits of the y-axis.}
  \item{userPrompt}{\bold{not used}}
  \item{save}{\bold{not used}}
  \item{xLeg}{x-coordinate for legend.}
  \item{yLeg}{y-coordinate for legend.}
  \item{yaxis.by}{increments along the y-axis to place tick marks.}
  \item{tcl.val}{tick length.}
  \item{B.col}{colour for spawning biomass.}
  \item{VB.col}{colour for vulnerable biomass.}
  \item{\dots}{\bold{not used}}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  AME: tried in separate file, but then changed that to
  lattice and wouldn't be good format for Arni's boxplots.

  Based on plotVBcatch (tweaking some). \cr
  \code{currentRes1} is local \code{currentRes}.
}
\seealso{
  \code{\link{plotVBcatch}}, \code{\link{plotBmcmcPOP}}
}
\keyword{hplot}


\name{plotVBcatch}
\alias{plotVBcatch}
\title{
  Plot Vulnerable Biomass
}
\description{
  Plot vulnerable biomass from posterior as boxplots, and
  add catch bars on same graph.
}
\usage{
plotVBcatch( obj, currentRes1=currentRes, 
             p=c(0.025,0.25,0.5,0.75,0.975),
             xyType="quantBox",
             lineType=c(3,2,1,2,3),
             refLines=NULL, xLim=NULL, yLim=NULL,
             userPrompt=FALSE, save=T, xLab = c(1939,1939),
             yLab=c(10000,220000),
             textLab=c("catch","vulnerable"),
             yaxis.by=10000, tcl.val=-0.2, ...)
}
\arguments{
  \item{obj}{MCMC data frame of \eqn{B} (\code{currentMCMC$B}).}
  \item{currentRes1}{list/scape object created by \code{\link{importCol2}}.}
  \item{p}{quantiles to use in \code{quantBox}.}
  \item{xyType}{type of plot (currently only uses \code{quantBox}).}
  \item{lineType}{line types to use in \code{quantBox}.}
  \item{refLines}{reference lines to add to plot.}
  \item{xLim}{limits of the x-axis.}
  \item{yLim}{limits of the y-axis.}
  \item{userPrompt}{\bold{not used}}
  \item{save}{\bold{not used}}
  \item{xLab}{x-coordinates for labels.}
  \item{yLab}{y-coordinates for labels.}
  \item{textLab}{text labels to display on plot.}
  \item{yaxis.by}{increments along the y-axis to place tick marks.}
  \item{tcl.val}{tick length.}
  \item{\dots}{additional arguments passed to the function \code{rect}.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  AME: This function is essentially a tweak of \code{plotBmcmcPOP}.
}
\seealso{
  \code{\link{plotBmcmcPOP}}, \code{\link{plotB2}}
}
\keyword{hplot}


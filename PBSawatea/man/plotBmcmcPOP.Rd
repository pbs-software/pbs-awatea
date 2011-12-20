\name{plotBmcmcPOP}
\alias{plotBmcmcPOP}
\title{
  Plot Spawning and Vulnerable Biomass
}
\description{
  Plot spawning and vulnerable biomass from posterior as boxplots, and
  add catch bars on same graph.
}
\usage{
plotBmcmcPOP( obj, currentRes1=currentRes, 
              p=c(0.025,0.25,0.5,0.75,0.975),
              xyType="quantBox",
              lineType=c(3,2,1,2,3),
              refLines=NULL, xLim=NULL, yLim=NULL,
              userPrompt=FALSE, save=TRUE, xLab=c(1939,1939,1939),
              yLab=c(10000,70000,170000),
              textLab=c("catch","spawning","vulnerable"),
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
   Combines ideas from \code{plt.quantBio} and \code{plotB2}. 
   Don't need lattice, just one figure, no panels.
   Vulnerable biomass has no posterior saved, which must be why 
   it's not been done before. Hmmm.... still worth seeing spawning though?

  Taking what is needed from \code{plt.quantBio}, this basically works: \cr
  \code{plt.quantBio(currentMCMC$B, xyType=rpType)}, \cr 
  though it creates 2x3 plots. The object should be the specific MCMC 
  posterior by year (so just a data frame), \emph{e.g.}, \code{currentMCMC$B}. \cr
  \code{currentRes1} is local \code{currentRes}.
}
\seealso{
  \code{\link{plotVBcatch}}, \code{\link{plotB2}}
}
\keyword{hplot}


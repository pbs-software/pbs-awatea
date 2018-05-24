\name{plotRmcmcPOP}
\alias{plotRmcmcPOP}
\title{
  Plot Recruitment Posterior Quantiles
}
\description{
  Plot recruitment posteriors quantiles as one graph over time.
}
\usage{
plotRmcmcPOP( obj, 
              p=get("quants5"),
              xyType="quantBox",
              lineType=c(3,2,1,2,3),
              refLines=NULL, xLim=NULL, yLim=NULL,
              userPrompt=FALSE, save=TRUE, tcl.val=-0.2,
              yaxis.by=10000, yLab="Recruitment", ...)
}
\arguments{
  \item{obj}{MCMC data frame of \eqn{R} (\code{currentMCMC$R}).}
  \item{p}{quantiles to use in \code{quantBox}.}
  \item{xyType}{type of plot (currently only uses \code{quantBox}).}
  \item{lineType}{line types to use in \code{quantBox}.}
  \item{refLines}{reference lines to add to plot.}
  \item{xLim}{limits of the x-axis.}
  \item{yLim}{limits of the y-axis.}
  \item{userPrompt}{\bold{not used}}
  \item{save}{\bold{not used}}
  \item{tcl.val}{tick length.}
  \item{yaxis.by}{increments along the y-axis to place tick marks.}
  \item{yLab}{label for the y-axis.}
  \item{\dots}{additional arguments passed to the function \code{rect}.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  AME: Plot recruitment posteriors quantiles as one graph over time. \cr
  Already have the full posterior densities done. \cr
  Using \code{plotBmcmcPOP} as template, but will be simpler as no extra stuff. 
  Probably not simplifying down as much as could due to time constraints. \cr
  Adding \code{yLab} and then using for exploitation plot also.
}
\seealso{
  \code{\link{plotB2}}
}
\keyword{hplot}


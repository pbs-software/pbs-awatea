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
   yaxis.by=10000, yLab="Recruitment", lang="e", ...)
}
\arguments{
  \item{obj}{\code{data.frame} -- MCMC data frame of \eqn{R} (\code{currentMCMC$R}).}
  \item{p}{\code{numeric} -- quantiles to use in \code{quantBox}.}
  \item{xyType}{\code{character} -- type of plot (currently only uses \code{quantBox}).}
  \item{lineType}{\code{numeric} -- line types to use in \code{quantBox}.}
  \item{refLines}{\code{numeric} -- reference lines to add to plot.}
  \item{xLim}{\code{numeric} -- limits of the x-axis.}
  \item{yLim}{\code{numeric} -- limits of the y-axis.}
  \item{userPrompt}{\code{logical} -- \bold{not used}}
  \item{save}{\code{logical} -- \bold{not used}}
  \item{tcl.val}{\code{numeric} -- tick length.}
  \item{yaxis.by}{\code{numeric} -- increments along the y-axis to place tick marks.}
  \item{yLab}{\code{character} -- label for the y-axis.}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
  \item{\dots}{additional arguments passed to the function \code{rect}.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-26}
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


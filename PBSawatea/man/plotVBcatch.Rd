\name{plotVBcatch}
\alias{plotVBcatch}
\title{
  Plot Vulnerable Biomass and Catch
}
\description{
  Plot vulnerable biomass from posterior as boxplots, and
  add catch bars on same graph.
}
\usage{
plotVBcatch(obj, currentRes1=currentRes,
   p=tcall(quants5),
   xyType="quantBox",
   lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, 
   xLab="Year",
   yLab="Catch and vulnerable biomass (t)",
   textLab=c("catch", "vulnerable"),
   yaxis.by=10000, tcl.val=-0.2,
   gear=1, lang="e", ...)
}
\arguments{
  \item{obj}{MCMC data frame of \emph{VB} (\code{currentMCMC$VB})}
  \item{currentRes1}{list/scape object created by \code{\link{importRes}}}
  \item{p}{quantiles to use in \code{quantBox}}
  \item{xyType}{type of plot (currently only uses \code{quantBox})}
  \item{lineType}{line types to use in \code{quantBox}}
  \item{refLines}{reference lines to add to plot}
  \item{xLim}{limits of the x-axis}
  \item{yLim}{limits of the y-axis}
  \item{userPrompt}{\bold{not used}}
  \item{save}{\bold{not used}}
  \item{xLab}{x-coordinates for labels}
  \item{yLab}{y-coordinates for labels}
  \item{textLab}{text labels to display on plot}
  \item{yaxis.by}{increments along the y-axis to place tick marks}
  \item{tcl.val}{tick length}
  \item{gear}{scalar: numeric gear code from subset 1 to Ngear}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
  \item{\dots}{additional arguments passed to the function \code{rect}}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-26}
}
\note{
  This function is essentially a tweak of \code{plotBmcmcPOP} (AME),
  which is currenly dysfunctional and not used (RH).
}
\seealso{
  \code{\link{plotB2}}
}
\keyword{hplot}


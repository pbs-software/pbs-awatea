\name{plotVBcatch}
\alias{plotVBcatch}
\title{
  Plot Vulnerable Biomass and Catch
}
\description{
  Plot the vulnerable biomass from posterior distributions as boxplots,
  and add catch bars in the same figure.
}
\usage{
plotVBcatch(obj, currentRes1=currentRes,
   p=tcall(quants5), xyType="quantBox", lineType=c(3,2,1,2,3),
   refLines=NULL, xLim=NULL, yLim=NULL,
   userPrompt=FALSE, save=TRUE, 
   xLab="Year",  yLab="Catch and vulnerable biomass (t)",
   textLab=c("catch", "vulnerable"),
   yaxis.by=10000, tcl.val=-0.1, cex.axis=1,
   gear=1, lang="e", ...)
}
\arguments{
  \item{obj}{\code{data.frame} -- MCMC data object of vulnerable biomass \emph{VB} over time (\code{currentMCMC$VB})}
  \item{currentRes1}{\code{list|scape} -- object created by \code{\link{importRes}}}
  \item{p}{\code{numeric} -- quantiles to use in \code{quantBox}}
  \item{xyType}{\code{character} -- type of plot (currently only uses \code{quantBox})}
  \item{lineType}{\code{numeric|character} -- line types to use in \code{quantBox}}
  \item{refLines}{\code{numeric} -- reference lines to add to plot}
  \item{xLim}{\code{numeric} -- limits of the x-axis}
  \item{yLim}{\code{numeric} -- limits of the y-axis}
  \item{userPrompt}{\code{logical} -- \bold{not used}}
  \item{save}{\code{logical} -- \bold{not used}}
  \item{xLab}{\code{numeric} -- x-coordinates for labels}
  \item{yLab}{\code{numeric} -- y-coordinates for labels}
  \item{textLab}{\code{character} -- text labels to display on plot}
  \item{yaxis.by}{\code{numeric} -- increments along the y-axis to place tick marks}
  \item{tcl.val}{\code{numeric} -- tick length}
  \item{cex.axis}{\code{numeric} -- size of axis labels}
  \item{gear}{\code{numeric} -- gear code from subset 1 to Ngear}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french)}
  \item{\dots}{\code{dots} -- additional arguments passed to the function \code{rect}}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Regional Headquarters (RHQ), Vancouver BC\cr
  Last modified \code{Rd: 2023-06-28}
}
\note{
  This function is essentially a tweak of \code{plotBmcmcPOP} (AME),
  which is currenly dysfunctional and not used (RH).
}
\seealso{
  \code{\link{plotB2}}
}
\keyword{hplot}


\name{plotAges}
\alias{plotAges}
\title{
  Plot MPD Fits to Proportions-at-Age
}
\description{
  Plot the MPD model fits from Awatea to age data (commercial or survey)
  using the scape function \code{plotCA}.
}
\usage{
plotAges(obj, what="c", maxcol=4, sexlab=c("Females","Males"), 
   ptypes=c("eps","png"), pngres=150, ...)
}
\arguments{
  \item{obj}{the MPD object created by the \code{importRes} function (e.g., \code{currentRes}).}
  \item{what}{character letter -- \code{"c"} for commercial proportions-at-age (pA) or \code{"s"} for survey pA.}
  \item{maxcol}{maximum number of columns for trellis plot}
  \item{sexlab}{label(s) to rename sex on the plot (e.g., a user may wish to rename "Unisex" to "Female").}
  \item{ptypes}{plot types for the figures: \code{eps} = encapsulated postscript, \code{png} = portable network graphics.}
  \item{pngres}{resolution (pixels/inch) for the \code{png} file.}
  \item{...}{additional arguments for \code{plotCA} and internal control.}
}
\details{
  The function drills down into the MPD object to extract data to pass to \pkg{scape}'s function \code{plotCA}.

  If all goes well, \code{.eps} and \code{.png} files will be produced for each series and sex 
  of the commercial \code{CAc} and/or survey \code{CAs} data in the MPD object.
}
\author{
  Rowan Haigh and Andrew M. Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link[scape]{plotCA}}, \code{\link[PBSawatea]{plt.bubbles}}, \code{\link[PBSawatea]{plt.recdev}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
% __ONLY ONE__ keyword per line

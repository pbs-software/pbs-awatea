\name{plotAges}
\alias{plotAges}
\title{
  Plot MPD Fits to Proportions-at-Age
}
\description{
  Plot the MPD model fits from Awatea to age data (commercial or survey)
  using the \pkg{scape} function \code{'plotCA'}.
}
\usage{
plotAges(obj, what="c", maxcol=5, sexlab=c("Females","Males"), 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"), ...)
}
\arguments{
  \item{obj}{\code{list} -- the MPD object created by the \code{'importRes'} function (e.g., \code{'currentRes'}).}
  \item{what}{\code{character} -- letter: \code{"c"} for commercial proportions-at-age (pA) or \code{"s"} for survey pA.}
  \item{maxcol}{\code{numeric} -- maximum number of columns for trellis plot.}
  \item{sexlab}{\code{character} -- label(s) to rename sex on the plot (e.g., a user may wish to rename \code{"Unisex"} to \code{"Female"}).}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{'eps'} = encapsulated postscript, \code{'png'} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{'png'} file.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
  \item{...}{additional arguments for \code{'plotCA'} and internal control.}
}
\details{
  The function drills down into the MPD object to extract data to pass
  to \pkg{scape}'s function \code{'plotCA'}.

  If desired, \code{.eps} and \code{.png} files will be produced
  for each series and sex of the commercial \code{'CAc'} and/or survey
  \code{'CAs'} data in the MPD object.
  However, the use of \code{'eps'} figures has largely been deprecated
  to build documents using MiKTeX's pdflatex engine, which requires
  \code{'png'} files.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2019-05-29}
}
\seealso{
In package \pkg{scape}:\cr
  \code{\link[scape]{plotCA}}\cr
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{plt.bubbles}},
  \code{\link[PBSawatea]{plt.recdev}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
% __ONLY ONE__ keyword per line

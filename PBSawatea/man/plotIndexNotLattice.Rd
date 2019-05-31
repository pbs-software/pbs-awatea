\name{plotIndexNotLattice}
\alias{plotIndexNotLattice}
\title{
  Plot Survey Indices without using Lattice
}
\description{
  Plot survey index series with error bars in a variety of ways (e.g., grouped, individually).
}
\usage{
plotIndexNotLattice(obj, main="", save=NULL,
   bar=1.96, ssnames=paste("Ser",1:9,sep=""),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"), ...)
}
\arguments{
  \item{obj}{\code{list} -- the MPD object created by the \code{'importRes'} function (default = \code{'currentRes'}).}
  \item{main}{\code{character} -- title for figure.}
  \item{save}{\code{deprecate} -- \emph{not used}}
  \item{bar}{\code{numeric} -- standard deviation of the normal distribution (1.96 is the approximate value of the 97.5 percentile point).}
  \item{ssnames}{\code{character} -- survey series names for figure labelling and image creation.}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{'eps'} = encapsulated postscript, \code{'png'} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{'png'} file.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
  \item{\dots}{\code{dots} -- \emph{not currently used}}
}
\value{
  Four series of portable network graphics (and encapsulated postscript, if needed) files: \cr
  \code{survIndSer.png......}{each survey panel focuses on the years of the survey;} \cr
  \code{survIndSer2.png.....}{each panel uses a fixed set of years that span all surveys;} \cr
  \code{survIndSer3.png.....}{one panel showing all series (incl. CPUE if used) normalised to their means;} \cr
  \code{survIndSer4-n.png...}{seperate files for each index series \emph{n}.}
}
\author{
First created by:\cr
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew M. Edwards}, PBS, DFO, Nanaimo BC

Modified by:\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2019-05-29}
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{plt.idx}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


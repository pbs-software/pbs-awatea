\name{plotCPUE}
\alias{plotCPUE}
\title{
  Plot CPUE with Error Bars
}
\description{
  Plot CPUE and fit with error bars.
}
\usage{
plotCPUE(obj, main="", save=NULL, bar=1.96, yLim=NULL, 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"), ...)
}
\arguments{
  \item{obj}{\code{data.frame} -- CPUE indices from Awatea's results file (\code{currentRes$CPUE}).}
  \item{main}{\code{character} -- title for figure.}
  \item{save}{\code{deprecate} -- \emph{not used}}
  \item{bar}{\code{numeric} -- standard deviation of the normal distribution (1.96 is the approximate value of the 97.5 percentile point).}
  \item{yLim}{\code{numeric} -- limits of the y-axis.}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{'eps'} = encapsulated postscript, \code{'png'} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{'png'} file.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
  \item{\dots}{\code{dots} -- \emph{not used}}
}
\value{
  A portable network graphics file: \cr
  \code{CPUEser.png} -- CPUE indices with error bars.
}
\author{
First created by:\cr
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew M. Edwards}, PBS, DFO, Nanaimo BC

Modified by:\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\note{
  Copied code from \code{plotIndexNotLattice}.
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{plotIndexNotLattice}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{plt.idx}
\alias{plt.idx}
\title{
  Plot Survey Index Residuals
}
\description{
  Plot the survey index residuals as a quantile-quantile plot 
  (see \code{\link{qqnorm}}).
}
\usage{
plt.idx(obj, main="Residuals", save=NULL, ssnames=paste("Ser",1:9,sep=""),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"), ...)
}
\arguments{
  \item{obj}{\code{data.frame} -- data with columns \code{Year}, \code{stdRes}, and \code{Fit}.}
  \item{main}{\code{character} -- title for the plot}
  \item{save}{\code{logical} -- if \code{TRUE}, save the figure to a raster file (\code{.png}).}
  \item{ssnames}{\code{character} -- survey series names for figure labelling and image creation.}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{"eps"} = encapsulated postscript, \code{"png"} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{.png} file.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
  \item{\dots}{\bold{not used}}
}
\details{
  QQ-plots show sample quantiles vs. theoretical quantiles.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-26}
}
\note{
  The \code{save} option has been disabled for some reason.
}
\seealso{
  \code{\link{plt.stdResids}}, \code{\link{plotIndexNotLattice}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


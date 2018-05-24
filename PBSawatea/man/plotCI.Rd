\name{plotCI}
\alias{plotCI}
\title{
  Plot Confidence Intervals
}
\description{
  Given a set of x and y values with upper and lower bounds,
  plot the points with the vertical limits.
}
\usage{
plotCI(x, y=NULL, ui, li, uiw=0.05, liw=uiw, clipNA=TRUE, gap=1, 
   col=par("col"), barcol=col, lwd=par("lwd"), lty=par("lty"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y}{\code{numeric} -- coordinates for the center of vertical limits. 
    If only \code{x} is provided, \code{x} becomes \code{y} and \code{x} defaults to \code{1:length(x)}.}
  \item{ui}{\code{numeric} -- vector of heights above the points \code{(x,y)}. Set to NULL or NA to omit upper bars.}
  \item{li}{\code{numeric} -- vector of heights below the points \code{(x,y)}. Set to NULL or NA to omit lower bars.}
  \item{uiw}{\code{numeric} -- width of the edge of the upper bar's head (inches).}
  \item{liw}{\code{numeric} -- width of the edge of the lower bar's foot (inches).}
  \item{clipNA}{\code{logical} -- if \code{TRUE}, remove \code{NA} values on either end of the points vector \code{(x,y)}.}
  \item{gap}{\code{numeric} -- vertical space between the points and the bars.}
  \item{col}{\code{character} -- colour of the points.}
  \item{barcol}{\code{character} -- colour of the bars.}
  \item{lwd}{\code{numeric} -- line width of bars.}
  \item{lty}{\code{numeric} -- line type for bars.}
  \item{\dots}{additional parameters for plotting (e.g., \code{bg} for background colour of points)}
}
\details{
  Used primarily to plot abundance indices and their confidence limits.
}
\value{
  Invisibly returns a list of \code{x} and \code{y}.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-05-23}
}
\note{
  This function offers a pale imitation of \code{gplots::plotCI}, and was created to reduce package dependenciesin \code{PBSawatea}.
}
\seealso{
\code{\link[PBSawatea]{plotACFs}},
\code{\link[PBSawatea]{plotCPUE}},
\code{\link[PBSawatea]{plotMeanAge}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")

\name{plt.stdResids}
\alias{plt.stdResids}
\title{
  Plot Diagnostics for Standardised Residuals
}
\description{
  Plot standardised residuals against year, fitted value, and theoretical residuals.
}
\usage{
plt.stdResids(obj, pct=c(5,25,50,75,95), 
     main=NULL, yLim=NULL, xLim=xLim, lang="e")
}
\arguments{
  \item{obj}{\code{data.frame} -- data with columns \code{Year}, \code{stdRes}, and \code{Fit}.}
  \item{pct}{\code{numeric} -- percentiles to display as horizontal lines on the quantile-quantile plot.}
  \item{main}{\code{character} -- title for the figure.}
  \item{yLim}{\code{numeric} -- limits of the y-axis.}
  \item{xLim}{\code{numeric} -- limits of the x-axis.}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
}
\details{
  Figure provides three panels of standardised residuals \emph{vs.} \cr
  (i) years, (ii) fitted or predicted values, and (iii) theoretical quantiles.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-26}
}
\seealso{
  \code{\link{plt.ssbVbCatch}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

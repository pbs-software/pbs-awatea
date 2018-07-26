\name{plotChains}
\alias{plotChains}
\title{Plot Cumulative Frequency of MCMC Chains}
\description{
  Plot cumulative fequency of \eqn{n} chains by partitioning one trace. \cr
  (Modified from the function \code{plotTracePOP}.)
}
\usage{
plotChains(mcmc, nchains=3, pdisc=0.1, axes=FALSE, same.limits=FALSE,
   between=list(x=axes,y=axes), div=1, span=1/4, log=FALSE,
   base=10, main=NULL, xlab=NULL, ylab=NULL, cex.main=1.2,
   cex.lab=1, cex.strip=0.8, cex.axis=0.8,
   las=0, tck=0.4, tick.number=5, lty.trace=1, lwd.trace=1,
   col.trace="grey", lty.median=1, lwd.median=1,
   col.median="black", lty.quant=2, lwd.quant=1,
   col.quant="black", plot=TRUE, probs=get("quants3"), lang="e", ...)
}
\arguments{
  \item{mcmc}{MCMC chain(s) as a vector, data frame or \code{mcmc} object.}
  \item{nchains}{number of chains to create from one trace.}
  \item{pdisc}{proportion of the initial trace to discard before creating chains.}
  \item{axes}{whether axis values should be plotted.}
  \item{same.limits}{whether panels should have same x-axis limits.}
  \item{between}{list with \code{x} and \code{y} indicating panel spacing.}
  \item{div}{denominator to shorten values on the y axis.}
  \item{span}{smoothness parameter (\bold{not used}).}
  \item{log}{whether values should be log-transformed.}
  \item{base}{logarithm base.}
  \item{main}{main title.}
  \item{xlab}{x-axis title.}
  \item{ylab}{y-axis title.}
  \item{cex.main}{size of main title.}
  \item{cex.lab}{size of axis labels.}
  \item{cex.strip}{size of strip labels.}
  \item{cex.axis}{size of tick labels.}
  \item{las}{orientation of tick labels: 0=parallel, 1=horizontal,
    2=perpendicular, 3=vertical.}
  \item{tck}{tick mark length.}
  \item{tick.number}{number of tick marks.}
  \item{lty.trace}{line type of trace.}
  \item{lwd.trace}{line width of trace.}
  \item{col.trace}{colour of trace.}
  \item{lty.median}{line type of median.}
  \item{lwd.median}{line width of median.}
  \item{col.median}{colour of median.}
  \item{lty.quant}{line type of quantile trace.}
  \item{lwd.quant}{line width of quantile trace.}
  \item{col.quant}{colour of quantile trace.}
  \item{plot}{whether to draw plot.}
  \item{probs}{quantile values for quantile trace.}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
  \item{\dots}{passed to \code{panel.trace} (\bold{not used}).}
}
\value{
  When \code{plot=TRUE}, a trellis plot is drawn and a data frame is
  returned, containing the data used for plotting. When
  \code{plot=FALSE}, a trellis object is returned.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-25}
}
\note{
  This idea stemmed from a discussion with PJS.
}
\seealso{
  \code{\link{plotTracePOP}}, \code{\link{plotDensPOP}}
}
% Graphics
\keyword{hplot}
% Statistics
\keyword{distribution}

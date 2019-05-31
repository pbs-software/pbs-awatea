\name{panelChains}
\alias{panelChains}
\title{Plot Cumulative Frequency of MCMC Split Chains}
\description{
  Plot cumulative frequency of \eqn{n} chains by partitioning one trace (set of samples). \cr
}
\usage{
panelChains(mcmc, nchains=3, pdisc=0.1, axes=FALSE, same.limits=FALSE,
   between=list(x=axes,y=axes), div=1, span=1/4, log=FALSE,
   base=10, main=NULL, xlab=NULL, ylab=NULL, cex.main=1.2,
   cex.lab=1, cex.strip=0.8, cex.axis=0.8,
   las=0, tck=0.4, tick.number=5, lty.trace=1, lwd.trace=1,
   col.trace="grey", lty.median=1, lwd.median=1,
   col.median="black", lty.quant=2, lwd.quant=1,
   col.quant="black", plot=TRUE, probs=tcall(quants3), lang="e", ...)
}
\arguments{
  \item{mcmc}{\code{matrix|data.frame} -- MCMC samples (by row) of parameters or derived quantities (by column).}
  \item{nchains}{\code{numeric} -- number of chains to create from one trace of MCMC samples.}
  \item{pdisc}{\code{numeric} -- proportion of the initial trace to discard before creating chains.}
  \item{axes}{\code{logical} -- if \code{TRUE}, plot axis values.}
  \item{same.limits}{\code{logical} -- if \code{TRUE}, panels should have same x-axis limits.}
  \item{between}{\code{deprecate} -- list with \code{x} and \code{y} indicating panel spacing (\emph{not used}).}
  \item{div}{\code{numeric} -- denominator to shorten values on the y axis.}
  \item{span}{\code{numeric} -- smoothness parameter (\bold{not used}).}
  \item{log}{\code{logical} -- if \code{TRUE}, log-transform values.}
  \item{base}{\code{numeric} -- logarithm base.}
  \item{main}{\code{character} -- main title.}
  \item{xlab}{\code{character} -- x-axis title.}
  \item{ylab}{\code{character} -- y-axis title.}
  \item{cex.main}{\code{numeric} -- size of main title.}
  \item{cex.lab}{\code{numeric} -- size of axis labels.}
  \item{cex.strip}{\code{numeric} -- size of strip labels.}
  \item{cex.axis}{\code{numeric} -- size of tick labels.}
  \item{las}{\code{numeric} -- orientation of tick labels: 0=parallel, 1=horizontal,
    2=perpendicular, 3=vertical.}
  \item{tck}{\code{numeric} -- tick mark length.}
  \item{tick.number}{\code{numeric} -- number of tick marks.}
  \item{lty.trace}{\code{numeric|character} -- line type of trace.}
  \item{lwd.trace}{\code{numeric} -- line width of trace.}
  \item{col.trace}{\code{numeric|character} -- colour of trace.}
  \item{lty.median}{\code{numeric|character} -- line type of median.}
  \item{lwd.median}{\code{numeric} -- line width of median.}
  \item{col.median}{\code{numeric|character} -- colour of median.}
  \item{lty.quant}{\code{numeric|character} -- line type of quantile trace.}
  \item{lwd.quant}{\code{numeric} -- line width of quantile trace.}
  \item{col.quant}{\code{numeric|character} -- colour of quantile trace.}
  \item{plot}{\code{logical} -- if \code{TRUE}, draw the plot.}
  \item{probs}{\code{numeric} -- quantile values for quantile trace.}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
  \item{\dots}{\code{dots} -- passed to \code{panel.trace} (\emph{not used}).}
}
\value{
  When \code{plot=TRUE}, a trellis plot is drawn and a data frame is
  returned, containing the data used for plotting. When
  \code{plot=FALSE}, a trellis object is returned.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2019-05-28}
}
\note{
  This idea stemmed from a discussion with PJS.
  Originally, the function was called \code{'plotChains'}.
}
\seealso{
  \code{\link{plotTracePOP}}, \code{\link{plotDensPOP}}
}
% Graphics
\keyword{hplot}
% Statistics
\keyword{distribution}

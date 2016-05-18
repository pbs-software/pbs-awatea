\name{plotTracePOP}
\alias{plotTracePOP}
\title{Plot MCMC Traces (AME Version)}
\description{
  Plot Markov-chain Monte Carlo traces. This is a diagnostic plot for
  deciding whether a chain shows unwanted trends by showing the trace 
  of moving quantiles. (Modified from the \pkg{plotMCMC} function
  \code{plotTrace}.)
}
\usage{
plotTracePOP(mcmc, axes=FALSE, same.limits=FALSE,
          between=list(x=axes,y=axes), div=1, span=1/4, log=FALSE,
          base=10, main=NULL, xlab=NULL, ylab=NULL, cex.main=1.2,
          cex.lab=1, cex.strip=0.8, cex.axis=0.8,
          las=0, tck=0.5, tick.number=5, lty.trace=1, lwd.trace=1,
          col.trace="grey", lty.median=1, lwd.median=1,
          col.median="black", lty.quant=2, lwd.quant=1,
          col.quant="black", plot=TRUE, probs=c(0.025, 0.5, 0.975),
          mpd=mcmc[1, ], ...)
}
\arguments{
  \item{mcmc}{MCMC chain(s) as a vector, data frame or \code{mcmc}
    object.}
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
  \item{mpd}{vector of MPD values (from which MCMC search starts).}
  \item{\dots}{passed to \code{panel.trace} (\bold{not used}).}
}
\value{
  When \code{plot=TRUE}, a trellis plot is drawn and a data frame is
  returned, containing the data used for plotting. When
  \code{plot=FALSE}, a trellis object is returned.
}
\note{
  The \code{Args} function from the \pkg{gdata} package is recommended
  for reviewing the arguments, instead of \code{args}.
}
\seealso{
  \code{\link[lattice]{xyplot}} and \code{\link[lattice]{panel.loess}}
  are the underlying drawing functions, and
  \code{\link[coda]{traceplot}} is a similar non-trellis plot.

  \code{\link[plotMCMC]{plotTrace}}, \code{\link[plotMCMC]{plotAuto}}, \code{\link[plotMCMC]{plotCumu}}, and
  \code{\link[plotMCMC]{plotSplom}} are diagnostic plots.

  \code{\link{plotDensPOP}}, \code{\link[plotMCMC]{plotDens}}, \code{\link[plotMCMC]{plotQuant}}, and 
  \code{\link{plotChains}} are posterior plots.

  \code{\link[plotMCMC]{plotMCMC-package}} gives an overview of the \pkg{plotMCMC} package.
}
% Graphics
\keyword{hplot}
% Statistics
\keyword{ts}

\name{plotDensPOP}
\alias{plotDensPOP}
\alias{plotDensPOPpars}
\alias{plotDensPOPparsPrior}
\title{Plot MCMC Density (AME Version)}
\description{
  Plot Markov-chain Monte Carlo density. This is an approximation of the
  posterior probability density function.
}
\usage{
plotDensPOP(mcmc, probs=c(0.025,0.975), points=FALSE, axes=TRUE,
     same.limits=FALSE, between=list(x=axes,y=axes), div=1,
     log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL,
     cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.7,
     las=0, tck=0.5, tick.number=5,
     lty.density=1, lwd.density=3, col.density="black",
     lty.median=2, lwd.median=1, col.median="darkgrey", 
     lty.outer=3, lwd.outer=1, col.outer="darkgrey", pch="|", 
     cex.points=1, col.points="black", plot=TRUE, MPD.height=0.04,...)

plotDensPOPpars(mcmc, probs=c(0.025,0.975), points=FALSE, 
     axes=TRUE, same.limits=FALSE, between=list(x=axes,y=axes), 
     div=1, log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
     cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.7, 
     las=0, tck=0.5, tick.number=5, 
     lty.density=1, lwd.density=3, col.density="black", 
     lty.median=2, lwd.median=1, col.median="darkgrey", 
     lty.outer=3, lwd.outer=1, col.outer="darkgrey", pch="|", 
     cex.points=1, col.points="black", plot=TRUE, MPD.height=0.04,...) 

plotDensPOPparsPrior(mcmc, probs=c(0.025,0.975), points=FALSE, 
     axes=TRUE, same.limits=FALSE, between=list(x=axes,y=axes), 
     div=1, log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
     cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.7, 
     las=0, tck=0.5, tick.number=5, 
     lty.density=1, lwd.density=3, col.density="black", 
     lty.median=2, lwd.median=1, col.median="darkgrey", 
     lty.outer=3, lwd.outer=1, col.outer="darkgrey", pch="|", 
     cex.points=1, col.points="black", plot=TRUE, MPD.height=0.04,...)
}
\arguments{
  \item{mcmc}{MCMC chain(s) as a vector, data frame or \code{mcmc} object.}
  \item{probs}{vector of outer quantiles to draw, besides the median.}
  \item{points}{whether data points should be plotted along the x axis.}
  \item{axes}{whether axis values should be plotted.}
  \item{same.limits}{whether panels should have same x-axis limits.}
  \item{between}{list with \code{x} and \code{y} indicating panel spacing.}
  \item{div}{denominator to shorten values on the x axis.}
  \item{log}{whether values should be log-transformed.}
  \item{base}{logarithm base.}
  \item{main}{main title.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{cex.main}{size of main title.}
  \item{cex.lab}{size of axis labels.}
  \item{cex.strip}{size of strip labels.}
  \item{cex.axis}{size of tick labels.}
  \item{las}{orientation of tick labels: 0=parallel, 1=horizontal, 2=perpendicular, 3=vertical.}
  \item{tck}{tick mark length.}
  \item{tick.number}{number of tick marks.}
  \item{lty.density}{line type of density curve.}
  \item{lwd.density}{line width of density curve.}
  \item{col.density}{colour of density curve.}
  \item{lty.median}{line type of median.}
  \item{lwd.median}{line width of median.}
  \item{col.median}{colour of median.}
  \item{lty.outer}{line type of outer quantiles.}
  \item{lwd.outer}{line width of outer quantiles.}
  \item{col.outer}{colour of outer quantiles.}
  \item{pch}{symbol for data points.}
  \item{cex.points}{size of data points.}
  \item{col.points}{colour of data points.}
  \item{plot}{whether to draw plot.}
  \item{MPD.height}{how far up to put MPD.}
  \item{\dots}{passed to \code{densityplot} and \code{panel.densityplot}.}
}
\details{
  The function \code{plotDensPOPpars} differs from \code{plotDensPOP} only 
  by a few tweaks to the internal list object \code{myscales}. \cr
  The function \code{plotDensPOPparsPrior} adds the priors automatically.
}
\value{
  When \code{plot=TRUE}, a trellis plot is drawn and a data frame is
  returned, containing the data used for plotting. When \code{plot=FALSE}, 
  a trellis object is returned.
	}
\note{
  The \code{Args} function from the \pkg{gdata} package is recommended
  for reviewing the arguments, instead of \code{args}.

  AME: edited \code{plotDens} function to have less whitesapce, 
  not repeat x-axis labels, and make y-axes the same scales. 
  Cannot just do it through the options.
  
  For Recruits and Biomass, use \code{plotDensPOPpars} for parameters.
  Tried y-axes the same scales, but 1973--1975 are so narrow that they
  make all the others really small: \code{same.limits=TRUE, ylim=c(0, 0.0005)}.
}
\seealso{
  \code{\link[lattice]{xyplot}} and
  \code{\link[lattice]{panel.densityplot}} are the underlying drawing
  functions, and \code{\link[coda]{densplot}} is a similar non-trellis
  plot.

  \code{\link[scapeMCMC]{plotTrace}}, \code{\link[scapeMCMC]{plotAuto}},
  \code{\link[scapeMCMC]{plotCumu}}, and \code{\link[scapeMCMC]{plotSplom}} are diagnostic
  plots.

  \code{\link[scapeMCMC]{plotDens}} and \code{\link[scapeMCMC]{plotQuant}} are posterior plots.

  \code{\link[scapeMCMC]{scapeMCMC-package}} gives an overview of the package.
}
% Graphics
\keyword{hplot}
% Statistics
\keyword{ts}

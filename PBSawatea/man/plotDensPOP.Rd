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
plotDensPOP(mcmc, probs=tcall(quants3)[c(1, 3)], points=FALSE, axes=TRUE,
     same.limits=FALSE, between=list(x=axes,y=axes), div=1,
     log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL,
     cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.7,
     las=0, tck=0.5, tick.number=5,
     lty.density=1, lwd.density=3, col.density="black",
     lty.median=2, lwd.median=1, col.median="darkgrey", 
     lty.outer=3, lwd.outer=1, col.outer="darkgrey", pch="|", 
     cex.points=1, col.points="black", plot=TRUE, MPD.height=0.04,
     mpd = mcmc[1, ], lang="e", ...)

plotDensPOPpars(mcmc, probs=tcall(quants3)[c(1, 3)], points=FALSE, 
     axes=TRUE, same.limits=FALSE, between=list(x=axes,y=axes), 
     div=1, log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
     cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.7, 
     las=0, tck=0.5, tick.number=5, 
     lty.density=1, lwd.density=3, col.density="black", 
     lty.median=2, lwd.median=1, col.median="darkgrey", 
     lty.outer=3, lwd.outer=1, col.outer="darkgrey", pch="|", 
     cex.points=1, col.points="black", plot=TRUE, MPD.height=0.04,
     mpd = mcmc[1, ], lang="e", ...) 

plotDensPOPparsPrior(mcmc, probs=tcall(quants3)[c(1, 3)], points=FALSE, 
     axes=TRUE, same.limits=FALSE, between=list(x=axes,y=axes), 
     div=1, log=FALSE, base=10, main=NULL, xlab=NULL, ylab=NULL, 
     cex.main=1.2, cex.lab=1, cex.strip=0.8, cex.axis=0.7, 
     las=0, tck=0.5, tick.number=5, 
     lty.density=1, lwd.density=3, col.density="black", 
     lty.median=2, lwd.median=1, col.median="darkgrey", 
     lty.outer=3, lwd.outer=1, col.outer="darkgrey", pch="|", 
     cex.points=1, col.points="black", plot=TRUE, MPD.height=0.04,
     mpd = mcmc[1, ], lang="e", ...)
}
\arguments{
  \item{mcmc}{\code{numeric} -- MCMC chain(s) as a vector, data frame or \code{mcmc} object.}
  \item{probs}{\code{numeric} -- vector of outer quantiles to draw, besides the median.}
  \item{points}{\code{logical} -- whether data points should be plotted along the x axis.}
  \item{axes}{\code{logical} -- whether axis values should be plotted.}
  \item{same.limits}{\code{logical} -- whether panels should have same x-axis limits.}
  \item{between}{\code{numeric} -- list with \code{x} and \code{y} indicating panel spacing.}
  \item{div}{\code{numeric} -- denominator to shorten values on the x axis.}
  \item{log}{\code{logical} -- whether values should be log-transformed.}
  \item{base}{\code{numeric} -- logarithm base}
  \item{main}{\code{character} -- main title}
  \item{xlab}{\code{character} -- x-axis label}
  \item{ylab}{\code{character} -- y-axis label}
  \item{cex.main}{\code{numeric} -- size of main title}
  \item{cex.lab}{\code{numeric} -- size of axis labels}
  \item{cex.strip}{\code{numeric} -- size of strip labels}
  \item{cex.axis}{\code{numeric} -- size of tick labels}
  \item{las}{\code{numeric} -- orientation of tick labels: 0=parallel, 1=horizontal, 2=perpendicular, 3=vertical}
  \item{tck}{\code{numeric} -- tick mark length}
  \item{tick.number}{\code{numeric} -- number of tick marks}
  \item{lty.density}{\code{numeric} -- line type of density curve}
  \item{lwd.density}{\code{numeric} -- line width of density curve}
  \item{col.density}{\code{character|numeric} -- colour of density curve}
  \item{lty.median}{\code{numeric} -- line type of median}
  \item{lwd.median}{\code{numeric} -- line width of median}
  \item{col.median}{\code{character|numeric} -- colour of median}
  \item{lty.outer}{\code{numeric} -- line type of outer quantiles}
  \item{lwd.outer}{\code{numeric} -- line width of outer quantiles}
  \item{col.outer}{\code{character|numeric} -- colour of outer quantiles}
  \item{pch}{\code{numeric} -- symbol for data points}
  \item{cex.points}{\code{numeric} -- size of data points}
  \item{col.points}{\code{character|numeric} -- colour of data points}
  \item{plot}{\code{logical} -- whether to draw plot}
  \item{MPD.height}{\code{numeric} -- how far up to put MPD}
  \item{mpd}{\code{numeric} -- vector of MPD values (from which MCMC search starts)}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
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
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-25}
}
\seealso{
  \code{\link[lattice]{xyplot}} and
  \code{\link[lattice]{panel.densityplot}} are the underlying drawing
  functions, and \code{\link[coda]{densplot}} is a similar non-trellis
  plot.

  \code{\link[plotMCMC]{plotTrace}}, \code{\link[plotMCMC]{plotAuto}},
  \code{\link[plotMCMC]{plotCumu}}, and \code{\link[plotMCMC]{plotSplom}} are diagnostic
  plots.

  \code{\link[plotMCMC]{plotDens}} and \code{\link[plotMCMC]{plotQuant}} are posterior plots.

  \code{\link[plotMCMC]{plotMCMC-package}} gives an overview of the package.
}
% Graphics
\keyword{hplot}
% Statistics
\keyword{ts}

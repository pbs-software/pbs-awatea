\name{panelTraces}
\alias{panelTraces}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot MCMC Traces in Panels
}
\description{
  Plot MCMC traces in panels for a given parameter, where each panel
  represents a separate run. That is, compare traces for a parameter
  among multiple model runs (e.g., base runs, sensitivity runs).
}
\usage{
panelTraces(mcmc, mpd=mcmc[1,], nchains=1, pdisc=0, axes=FALSE, 
   same.limits=FALSE, between=list(x=axes,y=axes), div=1, log=FALSE,
   base=10, main=NULL, xlab=NULL, ylab=NULL, cex.main=1.2, cex.lab=1.2,
   cex.strip=0.9, cex.axis=0.9, las=0, tck=0.4, tick.number=5,
   xfac=NULL, lang="e", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc}{\code{data.frame} -- MCMC samples (by row) of parameters (by column).}
  \item{mpd}{\code{numeric} -- vector of starting points for the MCMC simulations, usually the MPD.}
  \item{nchains}{\code{numeric} -- number of chains to create from one trace of MCMC samples.}
  \item{pdisc}{\code{numeric} -- proportion of the initial trace to discard before creating chains.}
  \item{axes}{\code{logical} -- if \code{TRUE}, plot axis values.}
  \item{same.limits}{\code{logical} -- if \code{TRUE}, panels should have same x-axis limits.}
  \item{between}{\code{deprecate} -- list with \code{x} and \code{y} indicating panel spacing (\emph{not used}).}
  \item{div}{\code{numeric} -- denominator to shorten values on the y axis.}
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
  \item{xfac}{\code{character} -- vector of labels for x-axis ticks (under quantile boxes).}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
  \item{\dots}{\code{dots} -- passed to sub function \code{'panel.trace'}.}
}
\details{
  Primarily used to compare how parameters compare across base runs or
  sensitivity runs when building the \sQuote{Model Results} appendix.
}
\value{
  Invisibly returns the data used by the function \code{'mochaLatte'} to plot the panels.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2019-05-29}
}
\note{
  Based on the function \code{'panelChains'}.
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{panelBoxes}},
  \code{\link[PBSawatea]{panelChains}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


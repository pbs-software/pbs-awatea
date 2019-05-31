\name{panelBoxes}
\alias{panelBoxes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Quantile Boxes in Panels
}
\description{
  Plot quantile boxes of MCMC samples for base or sensitivity runs in 
  panels that each describe a parameter and/or a derived quantity.
}
\usage{
panelBoxes(mcmc, nchains=9, pdisc=0, axes=FALSE, same.limits=FALSE,
   between=list(x=axes,y=axes), div=1, log=FALSE, base=10, main=NULL,
   xlab=NULL, ylab=NULL, xlim=c(0.25,9.75), 
   boxfill=paste0(rep(c("cyan","green","coral"),each=3),rep(1:3,3)),
   cex.main=1.2, cex.lab=1.2, cex.strip=0.9, cex.axis=0.9, las=0,
   tck=0.4, tick.number=5, xfac=paste0("B",1:nchains), outline=TRUE, 
   lang="e", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmc}{\code{data.frame} -- MCMC samples (by row) of parameters or derived quantities (by column).}
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
  \item{xlim}{\code{numeric} -- limits of the x-axis (mainly to add outer spacing to the outer boxes).}
  \item{boxfill}{\code{character} -- vector of colours to fill the quantile boxes.}
  \item{cex.main}{\code{numeric} -- size of main title.}
  \item{cex.lab}{\code{numeric} -- size of axis labels.}
  \item{cex.strip}{\code{numeric} -- size of strip labels.}
  \item{cex.axis}{\code{numeric} -- size of tick labels.}
  \item{las}{\code{numeric} -- orientation of tick labels: 0=parallel, 1=horizontal,
    2=perpendicular, 3=vertical.}
  \item{tck}{\code{numeric} -- tick mark length.}
  \item{tick.number}{\code{numeric} -- number of tick marks.}
  \item{xfac}{\code{character} -- vector of labels for x-axis ticks (under quantile boxes).}
  \item{outline}{\code{logical} -- if \code{TRUE}, add outliers to the quantile plots.}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
  \item{\dots}{\code{dots} -- passed to subfunction \code{'panel.box'}.}
}
\details{
  Primarily used to compare how parameters or derived quantities 
  compare across base runs or sensitivity runs when building the 
  \sQuote{Model Results} appendix.
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
  \code{\link[PBSawatea]{panelChains}},
  \code{\link[PBSawatea]{panelTraces}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


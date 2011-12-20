\name{plotB2}
\alias{plotB2}
\title{Plot Biomass, Recruitment, and Landings (AME Version)}
\description{
  Plot \code{scape} model predicted biomass, stock recruitment, and
  landings.

  AME: This is an alteration of Arni Magnussons \code{plotB} function to accommodate
  PJS's request not to show biomass prior to fishery and survey indices period.
}
\usage{
plotB2(model, what="d", series=NULL, years=NULL, axes=TRUE, div=1,
      legend="bottom", main="", xlab="", ylab="", cex.main=1.2,
      cex.legend=1, cex.lab=1, cex.axis=0.8, las=1,
      tck=c(1,what=="d")/2, tick.number=5, lty.grid=3, col.grid="white",
      pch=16, cex.points=0.8, col.points="black", lty.lines=1:3,
      lwd.lines=2, col.lines="black", ratio.bars=3, col.bars="grey",
      plot=TRUE, \dots)
}
\arguments{
  \item{model}{fitted \code{scape} model.}
  \item{what}{what to plot: \code{"d"}[efault], \code{"s"}[tock
    recruitment], or \code{"l"}[andings].}
  \item{series}{vector of strings indicating which column names in
    \code{model$B} data frame to plot (all by default).}
  \item{years}{vector of numbers indicating which years to include (all
    by default).}
  \item{axes}{whether to plot axis values.}
  \item{div}{denominator to shorten values on the y axis, or a vector
    with two elements referring to x and y axis.}
  \item{legend}{legend location: \code{"bottom"}, \code{"left"},
    \code{"top"}, \code{"right"}, or \code{""} to suppress legend.}
  \item{main}{main title.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{cex.main}{size of main title.}
  \item{cex.legend}{size of legend text.}
  \item{cex.lab}{size of axis labels.}
  \item{cex.axis}{size of tick labels.}
  \item{las}{orientation of tick labels: 0=parallel, 1=horizontal,
    2=perpendicular, 3=vertical.}
  \item{tck}{tick mark length.}
  \item{tick.number}{number of tick marks.}
  \item{lty.grid}{line type of gridlines.}
  \item{col.grid}{color of gridlines.}
  \item{pch}{symbol for points.}
  \item{cex.points}{size of points.}
  \item{col.points}{color of points.}
  \item{lty.lines}{line type of main lines.}
  \item{lwd.lines}{line width of main lines.}
  \item{col.lines}{color of main lines.}
  \item{ratio.bars}{width of bars.}
  \item{col.bars}{color of bars.}
  \item{plot}{whether to draw plot.}
  \item{\dots}{passed to \code{xyplot} and \code{panel.superpose}.}
}
\details{
  The \code{"d"}[efault] plot shows spawning biomass and vulnerable
  biomass as lines, and landings as bars, on the same scale.
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
  \code{\link[lattice]{xyplot}}, \code{\link[lattice]{panel.barchart}},
  and \code{\link[lattice]{panel.superpose}} are the underlying drawing
  functions.

  \code{\link[scape]{plotCA}}, \code{\link[scape]{plotCL}},
  \code{\link[scape]{plotIndex}}, \code{\link[PBSawatea]{plotIndex2}}
  and \code{\link[scape]{plotLA}} plot model fit and data.

  \code{\link[scape]{plotB}}, \code{\link[scape]{plotN}}, 
  and \code{\link[scape]{plotSel}} plot derived quantities.

  \code{\link[scape]{scape-package}} gives an overview of the \pkg{scape} package.
}
% Graphics
\keyword{hplot}


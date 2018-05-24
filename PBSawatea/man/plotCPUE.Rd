\name{plotCPUE}
\alias{plotCPUE}
\title{
  Plot CPUE and Add Error Bars
}
\description{
  Plot CPUE and fit with error bars.
}
\usage{
plotCPUE(obj, main="", save=NULL, bar=1.96, yLim=NULL, 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, ...)
}
\arguments{
  \item{obj}{data frame of CPUE indices from Awatea's results file
    (\emph{e.g.}, \code{currentRes$CPUE}).}
  \item{main}{\code{character} -- title for figure}
  \item{save}{\bold{not used}}
  \item{bar}{\code{numeric} -- standard deviation of the normal distribution 
    (1.96 is the approximate value of the 97.5 percentile point).}
  \item{yLim}{\code{numeric} -- limits of the y-axis.}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{eps} = encapsulated postscript, 
     \code{png} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{png} file.}
  \item{\dots}{\bold{not used}}
}
\value{
  A postscript file: \cr
  \code{CPUEser.eps...}{CPUE indices with error bars.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Copied code from \code{plotIndexNotLattice}.
}
\seealso{
  \code{\link{plotIndexNotLattice}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

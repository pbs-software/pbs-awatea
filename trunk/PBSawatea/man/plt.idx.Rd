\name{plt.idx}
\alias{plt.idx}
\title{
  Plot Survey Index Residuals
}
\description{
  Plot the survey index residuals as a quantile-quantile plot 
  (see \code{\link{qqnorm}}).
}
\usage{
plt.idx(obj, main="Residuals", save=NULL, ssnames=paste("Ser",1:9,sep=""), ...)
}
\arguments{
  \item{obj}{a data frame with columns \code{Year}, \code{stdRes}, and \code{Fit}.}
  \item{main}{title for the plot.}
  \item{save}{logical: if \code{TRUE}, save the figure to a raster file (\code{.png}).}
  \item{ssnames}{survey series names for figure labelling and image creation.}
  \item{\dots}{\bold{not used}}
}
\details{
  QQ-plots show sample quantiles vs. theoretical quantiles.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  The \code{save} option has been disabled for some reason.
}
\seealso{
  \code{\link{plt.stdResids}}, \code{\link{plotIndexNotLattice}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


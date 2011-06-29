\name{panLegend}
\alias{panLegend}
\title{
  Place a Legend in a Figure using Relative (0:1) Coordinates
}
\description{
  Place a legend in a figure by first setting the coordinate space 
  to lie between 0 and 1: \cr \code{par(usr=c(0,1,0,1))}.
}
\usage{
panLegend(x, y, legTxt, ...)
}
\arguments{
  \item{x}{relative x-coordinate.}
  \item{y}{relative y-coordinate.}
  \item{legTxt}{legend text to add to figure.}
  \item{\dots}{additional arguments sent to function \code{legend}.}
}
\note{
  Currently, this function does not reset the coordinate space to the original. \cr
  Use \code{\link{addLegend}} instead.
}
\seealso{
  \code{\link{addLabel}}, \code{\link{addLegend}}
}
\keyword{utilities}


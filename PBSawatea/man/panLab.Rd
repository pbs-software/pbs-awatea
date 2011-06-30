\name{panLab}
\alias{panLab}
\title{
  Write Text to Figure in Relative (0:1) Coordinates
}
\description{
  Write text to a figure by first setting the coordinate space 
  to lie between 0 and 1: \cr \code{par(usr=c(0,1,0,1))}.
}
\usage{
panLab(x, y, txt, ...)
}
\arguments{
  \item{x}{relative x-coordinate.}
  \item{y}{relative y-coordinate.}
  \item{txt}{text to add to figure.}
  \item{\dots}{additional arguments sent to function \code{text}.}
}
\note{
  Currently, this function does not reset the coordinate space to the original. \cr
  Use \code{\link[PBSmodelling]{addLabel}} instead.
}
\seealso{
  \code{\link[PBSmodelling]{addLabel}}, \code{\link[PBSmodelling]{addLegend}}
}
\keyword{utilities}


\name{allEqual}
\alias{allEqual}
\title{
  Are All Values Equal to the First?
}
\description{
  A short-cut function for \code{all(x==x[1])}, which asks are all values
  in \eqn{x} equl to the first value in \eqn{x}.
}
\usage{
allEqual(x)
}
\arguments{
  \item{x}{vector of values.}
}
\value{
\code{TRUE} or \code{FALSE}
}
\seealso{
  \code{\link[base]{all}}, \code{\link[PBSmodelling]{clearAll}}, \code{\link[PBSmodelling]{clipVector}}
}
\keyword{logic}


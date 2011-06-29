\name{getYrIdx}
\alias{getYrIdx}
\title{
  Select a Subset of Years for Plotting
}
\description{
  Select a subset of years for which many years are available. 
  The default is to select 5-year increments.
}
\usage{
getYrIdx(yrNames, mod=5)
}
\arguments{
  \item{yrNames}{vector (character or numeric) of years.}
  \item{mod}{select the years modulo \code{mod}.}
}
\value{
  Subset of input years that are modulo \code{mod}.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link{findPat}}, \code{\link{pad0}}
}
\keyword{utilities}


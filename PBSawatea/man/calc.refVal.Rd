\name{calc.refVal}
\alias{calc.refVal}
\title{
  Calculate Reference Value for Performance Measure
}
\description{
  Calculate the reference value for performance measures.
}
\usage{
calc.refVal(obj, refYrs, fun=mean)
}
\arguments{
  \item{obj}{scape Biomass matrix with \eqn{n} rows and \eqn{m} columns, \cr
    where \eqn{n} = number of MCMC samples, and \eqn{m} = number of years.}
  \item{refYrs}{numeric years in reference period.}
  \item{fun}{the function to apply to reference period \eqn{i}.}
}
\value{
  Returns a vector of length \code{nrow(obj)} reference values.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link{calc.projExpect}}, \code{\link{findTarget}}
}
\keyword{manip}


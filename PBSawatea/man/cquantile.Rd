\name{cquantile}
\alias{cquantile}
\alias{cquantile.vec}
\title{
  Running Quantile
}
\description{
  Creates a set of running quantiles from MCMC traces. \cr
  (Uses subfunction found in \pkg{coda}'s function \code{cumuplot}.)
}
\usage{
cquantile(z, probs)
cquantile.vec(z, prob)
}
\arguments{
  \item{z}{an MCMC object.}
  \item{probs}{vector of quantiles.}
  \item{prob}{single quantile.}
}
\value{
  \code{cquantile.......}{running quantile matrix} \cr
  \code{cquantile.vec...}{running quantile vector}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Arni Magnusson describes a running quatile as: \cr
  \dQuote{\emph{the evolution of the sample quantiles as a function of the number of iterations}}
}
\seealso{
  \code{\link[coda]{cumuplot}}, \code{\link{plotTracePOP}}
}
\keyword{arith}
\keyword{ts}% __ONLY ONE__ keyword per line


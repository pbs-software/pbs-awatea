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
  \item{z}{\code{numeric} -- an MCMC object}
  \item{probs}{\code{numeric} -- vector of quantiles}
  \item{prob}{\code{numeric} -- single quantile}
}
\value{
  \tabular{ll}{
    \code{cquantile}     \tab running quantile matrix\cr
    \code{cquantile.vec} \tab running quantile vector
  }
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC\cr
  Last modified \code{Rd: 2025-06-17}
}
\note{
  Arni Magnusson describes a running quantile as: \cr
  \dQuote{\emph{the evolution of the sample quantiles as a function of the number of iterations}}
}
\seealso{
In package \pkg{coda}:\cr
  \code{\link[coda]{cumuplot}}\cr
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{plotTracePOP}}
}
\keyword{arith}
\keyword{ts}% __ONLY ONE__ keyword per line
\concept{utilFuns.r}


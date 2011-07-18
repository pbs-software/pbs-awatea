\name{tabSAR}
\alias{tabSAR}
\title{
  Make Files of Reference Point Tables
}
\description{
  Generate comma-delimited, two-dimensional output tables (\code{.csv})
  from objects reporting probabilities and times for DFO reference points 
  and COSEWIC critera.
}
\usage{
tabSAR(models=paste("input-ymr",pad0(c(29,30),2),pad0(1,2),sep="."),
    pnam="refProbs3Gen90", tnam=c("Ttab0.5","Ttab0.8","Ttab0.95"), 
    cats=seq(0, 2500, 500), digits=2)
}
\arguments{
  \item{models}{names of binary system files that store the decision tables.}
  \item{pnam}{name of list object containing matrices of reference probabilities.}
  \item{tnam}{names of matrices reporting times to reach reference points/criteria.}
  \item{cats}{catch strategies (subset) to report in output tables.}
  \item{digits}{number of digits to retain after the decimal.}
}
\details{
  At the moment, the reference point objects are a list of matrices for 
  the probabilities and individual matrices for the times to reach the 
  references. The latter should also be a list of matrices to parallel the 
  probability object.
}
\value{
  No value is returned in R; however, comma-delimited text files (\code{.csv})
  are created in the user's working directory. The probaility tables have the suffix 
  \code{_prob.csv} while the time-to-target tables have the suffix \code{_targ.csv}.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanimo BC
}
\note{
  The subset of the master decision tables in the stock assessment are primarily used
  in the SAR (Science Advisory Report) for managers.
}
\seealso{
  \code{\link[PBSawatea]{refPoints}}, \code{\link[PBSawatea]{findTarget}}, \cr
  \code{\link[base]{cat}}, \code{\link[base]{round}}, 
  \code{\link[PBSmodelling]{pad0}}, \code{\link[PBSmodelling]{show0}}
}
\keyword{IO}
\keyword{file}% __ONLY ONE__ keyword per line


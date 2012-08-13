\name{importEva}
\alias{importEva}
\title{
  Import Awatea Hessian Eigenvalues
}
\description{
  Import a vector of eigenvalues of the Hessian from an MPD analysis file (\code{Awatea.eva}).
}
\usage{
importEva(eva.file)
}
\arguments{
  \item{eva.file}{vector of hessian eigenvalues; first element should be a positive definite number.}
}
\details{
   The eigenvalues of the Hessian represent the 2nd derivatives of the negative log-likelihood function. 
  If they are all positive it indicates a minimum.
  The first element of this vector should be positive definite before an MCMC is run.
}
\value{
  The output is a list object with only one component:
  \item{eva}{a numeric vector representing the eigenvalues of the Hessian.}
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link[PBSawatea]{importCor}}, \code{\link[PBSawatea]{importLik}}, \code{\link[PBSawatea]{importPar}},
  \code{\link[PBSawatea]{importRes}}, \code{\link[PBSawatea]{importStd}}
}
\keyword{file}
\keyword{interface}


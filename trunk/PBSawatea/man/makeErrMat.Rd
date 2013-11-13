\name{makeErrMat}
\alias{makeErrMat}
\title{
  Make Ageing Error Matrix for Awatea
}
\description{
  Make a simple ageing error matrix for Awatea.
}
\usage{
makeErrMat(N=60, ondiag=0.8, offdiag=0.1, corner=0.9)
}
\arguments{
  \item{N}{numeric scalar indicating number of age classes, which determines the dimension of the matrix. }
  \item{ondiag}{numeric value to appear along the matrix diagonal. }
  \item{offdiag}{numeric value to appear one cell to the left and right of the matrix diagonal. }
  \item{corner}{numeric value to appear in the top left and bottom right corners of the matrix. }
}
\value{
  Simple symmetric ageing error matrix.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link[PBStools]{plotProp}}, \code{\link[PBStools]{weightBio}}
}
\keyword{array}


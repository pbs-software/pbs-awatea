\name{makeCmat}
\alias{makeCmat}
\alias{makeRmat}
\title{
  Make Column (Row) Matrix Using a Vector
}
\description{
  Covert a vector into a one-column (or one-row) matrix.
}
\usage{
makeCmat(x, colname = "Y")
makeRmat(x, rowname = "Y")
}
\arguments{
  \item{x}{vector (numeric or character) to convert to a matrix with dimensions (\emph{n},1) or (1,\emph{n}).}
  \item{colname}{name of the single column in the matrix}
  \item{rowname}{name of the single row in the matrix}
}
\details{
  R has a tendency to lose dimensionality when a matrix subset has one column or one row.
  The dimensionality can be retained using the argument \code{drop=FALSE};
  however, there are time when a user cannot control operations using this trick (e.g., reading in a file).
  The functions \code{makeCmat} and \code{makeRmat} can restore dimensionality to a vector when needed.
}
\value{
  A matrix with one column or one row, depending on the function.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link[PBSawatea]{runADMB}}
}
\keyword{manip}

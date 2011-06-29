\name{stdRes.index}
\alias{stdRes.index}
\title{
  Calculate Standardised Residuals for Abundance Indices
}
\description{
  Calculate the standardised residuals for commercial and survey indices.
}
\usage{
stdRes.index(obj, label=NULL, prt=TRUE)
}
\arguments{
  \item{obj}{data frame of observed and fitted index values from Awatea's results file (\code{.res}).}
  \item{label}{general label for the output.}
  \item{prt}{logical: if \code{TRUE}, print the results.}
}
\value{
  Input data frame with additional column of standardised residuals.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link{stdRes.CA}}, \code{\link{importCol2}}, \code{\link{reweight}}
}
\keyword{distribution}


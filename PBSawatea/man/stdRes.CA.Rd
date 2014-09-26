\name{stdRes.CA}
\alias{stdRes.CA}
\title{
  Calculate Standardised Residuals for Robust Normal Likelihood
}
\description{
  Calculate the standardised residuals for Awatea's implementation 
  of the Fournier \emph{robustified} normal likelihood for 
  proportions-at-length.\cr
  Based on PJS's summary of the CASAL document and ACH's change to length.
}
\usage{
stdRes.CA(obj, trunc=3, myLab="Age Residuals", prt=TRUE)
}
\arguments{
  \item{obj}{scape/list object of Awatea's results file (\code{.res}).}
  \item{trunc}{maximum standardised residual; values greater than this are set to \code{trunc}.}
  \item{myLab}{general label for the output.}
  \item{prt}{logical: if \code{TRUE}, print the results.}
}
\value{
  List object of standardised residuals.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link{stdRes.index}}, \code{\link{importRes}}, \code{\link{reweight}}
}
\keyword{robust}
\keyword{distribution}


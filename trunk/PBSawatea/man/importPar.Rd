\name{importPar}
\alias{importPar}
\title{
  Import Awatea Parameter File
}
\description{
  Import an Awatea file of all parameters resulting from an MPD minimisation.
}
\usage{
importPar(par.file)
}
\arguments{
  \item{par.file}{parameter file (e.g., \code{Awatea.par}).}
}
\details{
  Extracts the data from character vectors and makes various data objects.
}
\value{
  The output is a list object with the following components:
  \item{par}{character vector representing the file line-by-line}
  \item{npars}{number of parameters reported in header}
  \item{fval}{objective function value reported in header}
  \item{maxgrad}{maximum gradient component reported in header}
  \item{parameters}{numerous reported parameters (e.g., \code{R0}) that will depend on the model}
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link[PBSawatea]{importCor}}, \code{\link[PBSawatea]{importRes}}, \code{\link[PBSawatea]{importStd}}
}
\keyword{file}
\keyword{interface}


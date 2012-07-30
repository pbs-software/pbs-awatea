\name{importCor}
\alias{importCor}
\title{
  Import Awatea Correlation File
}
\description{
  Import an Awatea correlation file that results from an MPD minimisation.
}
\usage{
importCor(cor.file)
}
\arguments{
  \item{cor.file}{correlation file (e.g., \code{Awatea.cor}).}
}
\details{
  Extracts the data from character vectors and makes various data objects.
}
\value{
  The output is a list object with the following components:
  \item{cfile}{ character vector representing the file line-by-line}
  \item{cor}{data.frame of the correlation file}
  \item{cor.mat}{matrix representing the correlation matrix only}
  \item{index}{character vector used in the row and column names of \code{cor} and \code{cor.mat}}
  \item{cor.name}{Awatea parameter names in the correlation file}
  \item{cor.value}{Awatea parameter values in the correlation file}
  \item{cor.std.dev}{Awatea parameter standard deviations in the correlation file}
  \item{hessian_log_determinant}{log of the determinant of the hessian reported in the header}
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link[PBSawatea]{importPar}}, \code{\link[PBSawatea]{importRes}}, \code{\link[PBSawatea]{importStd}}
}
\keyword{file}
\keyword{interface}


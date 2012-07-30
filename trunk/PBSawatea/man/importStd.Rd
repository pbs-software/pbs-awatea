\name{importStd}
\alias{importStd}
\title{
  Import Awatea Estimated Parameter File
}
\description{
  Import an Awatea file of estimated parameters resulting from an MPD minimisation.
}
\usage{
importStd(std.file, vnam="name")
}
\arguments{
  \item{std.file}{parameter file (e.g., \code{Awatea.std}).}
  \item{vnam}{field name listing parameter names in data file.}
}
\details{
  Extracts the data from character vectors and makes various data objects.
}
\value{
  The output is a list object with the following components:
  \item{std}{data.frame representation of the \code{.std} file}
  \item{parameters}{data.frame for each of the estimated parameters in \code{std}}
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link[PBSawatea]{importCor}}, \code{\link[PBSawatea]{importPar}}, \code{\link[PBSawatea]{importRes}}
}
\keyword{file}
\keyword{interface}


\name{load.allResFiles}
\alias{load.allResFiles}
\title{
  Load All Awatea .res Files
}
\description{
  Load all Awatea \code{.res} files in the working directory into a list object.
}
\usage{
load.allResFiles(resList = NULL)
}
\arguments{
  \item{resList}{AME: sets directory to path above current and 
    sets the pattern to \code{"results.dat$"}; probably deprecated.}
}
\value{
  List of multiple calls to \code{\link{importCol2}}.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  If deprecated, remove from package \pkg{PBSawatea}.
}
\seealso{
  \code{\link{importCol2}}
}
\keyword{list}


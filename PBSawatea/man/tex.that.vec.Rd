\name{tex.that.vec}
\alias{tex.that.vec}
\title{
  Tex That Vector
}
\description{
  Convert a vector to a phrase \sQuote{x, y, and z} for use in Sweave.
}
\usage{
tex.that.vec(vec, simplify=TRUE)
}
\arguments{
  \item{vec}{\code{numeric|character} -- an R vector to summarise in Sweave.}
  \item{simplify}{\code{logical} -- if \code{TRUE}, remove duplicated values and sort the unique values.}
}
\details{
  Primarily a utility function for use in Sweave where one wishes to convert a 
  vector of items into a phrase for a beautiful sentence.
}
\value{
	Returns a single-element character string that summarises the vector.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-05-24}
}
\note{
  This function was a pre-cursor to the more advanced \code{\link[PBStools]{texThatVec}}
  in \code{PBStools}; however, \code{PBSawatea} currently only needs \code{PBSmodelling}
  in an effort to reduce dependencies.
  
  \emph{Aside}: Dependence only on \code{PBSmodelling} is also the reason that the
  \code{PBStools} function \code{.flush.cat} is reproduced in \code{PBSawatea} as
  \code{.flash.cat}.
}
\seealso{
  \code{\link[PBSawatea]{panLab}},
  \code{\link[PBSawatea]{readAD}},
  \code{\link[PBSawatea]{tabSAR}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")

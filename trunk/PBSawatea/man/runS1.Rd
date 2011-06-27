\name{runS1}
\alias{runS1}
\title{
  Wrapper to function <runSweave>
}
\description{
  This small utility function simply provides a wrapper to the function 
  \code{runSweave} so that multiple documents can be built at once.
}
\usage{
runS1(runs=1, rewts=0:6, cpue=FALSE)
}
\arguments{
  \item{runs}{the run number(s).}
  \item{rewts}{the reweight number(s).}
  \item{cpue}{logical: if \code{TRUE}, retain the CPUE figures in the Sweave file.}
}
\details{
	Essentially loops through \code{\link{runSweave}} using \code{(i in runs)} and \code{(j in rewts)}.
}
\value{
  Produces multiple Sweave documents and PDF files.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  A fairly recent Seave file called \code{ymrrun-master.Snw} can be found in 
  the library directory: \cr
  \code{.../R/.../library/PBSawatea/snw}
}
\seealso{
  \code{\link{runSweave}}, \code{\link{runADMB}}
}
\keyword{character}
\keyword{manip}


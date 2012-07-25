\name{runMPD}
\alias{runMPD}
\title{
  Wrapper to Function <runSweave>
}
\description{
  This small utility function simply provides a wrapper to the function 
  \code{runSweave} so that multiple documents can be built at once.
}
\usage{
runMPD(prefix=c("spp","area"), runs=1, rwts=0, ...)
}
\arguments{
  \item{prefix}{character vector of \code{filename} prefix components delimited by \code{delim}.}
  \item{runs}{the run number(s).}
  \item{rwts}{the reweight number(s).}
  \item{...}{arguments passed on to the function \code{runSweave}.}
}
\details{
	Essentially loops through \code{\link{runSweave}} using \code{(i in runs)} and \code{(j in rwts)}.
}
\value{
  Produces multiple Sweave documents and PDF files.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  A recent Seave file called \code{ymrrun-master.Snw} can be found in 
  the library directory: \cr
  \code{.../R/.../library/PBSawatea/snw}
}
\seealso{
  \code{\link{runSweave}}, \code{\link{runSweaveMCMC}}, \code{\link{runMCMC}}
}
\keyword{character}
\keyword{manip}


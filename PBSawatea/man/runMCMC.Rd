\name{runMCMC}
\alias{runMCMC}
\title{
  Wrapper to Function <runSweaveMCMC>
}
\description{
  This small utility function simply provides a wrapper to the function 
  \code{runSweaveMCMC} so that multiple documents can be built at once.
}
\usage{
runMCMC(prefix=c("spp","area"), runs=1, rwts=0, ...)
}
\arguments{
  \item{prefix}{\emph{character} vector of \code{filename} prefix components delimited by \code{delim}.}
  \item{runs}{\emph{numeric} vector of run number(s).}
  \item{rwts}{\emph{numeric} vector of reweight number(s).}
  \item{...}{arguments passed on to the function \code{runSweaveMCMC}.}
}
\details{
	Essentially loops through \code{\link{runSweaveMCMC}} using \code{(i in runs)} and \code{(j in rewts)}.
}
\value{
  Produces multiple Sweave documents and PDF files for MCMCs.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link{runSweaveMCMC}}, \code{\link{runSweave}}, \code{\link{runMPD}}
}
\keyword{character}
\keyword{manip}


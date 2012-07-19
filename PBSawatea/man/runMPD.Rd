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
runMPD(strSpp="XYZ", prefix=c("spp","area"), 
    runs=1, rewts=0:6, cpue=FALSE, delim="-")
}
\arguments{
  \item{strSpp}{three-letter code that identifies the species.}
  \item{prefix}{character vector of \code{filename} prefix components delimited by \code{delim}.}
  \item{runs}{the run number(s).}
  \item{rewts}{the reweight number(s).}
  \item{cpue}{logical: if \code{TRUE}, retain the CPUE figures in the Sweave file.}
  \item{delim}{character that delimits the components of the \code{filename} prefix.}
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
  A recent Seave file called \code{ymrrun-master.Snw} can be found in 
  the library directory: \cr
  \code{.../R/.../library/PBSawatea/snw}
}
\seealso{
  \code{\link{runSweave}}, \code{\link{runSweaveMCMC}}, \code{\link{runMCMC}}
}
\keyword{character}
\keyword{manip}


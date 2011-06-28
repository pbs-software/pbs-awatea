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
runMCMC(runs=7, rewts=0:6, cpue=FALSE, estM=TRUE)
}
\arguments{
  \item{runs}{the run number(s).}
  \item{rewts}{the reweight number(s).}
  \item{cpue}{logical: if \code{TRUE}, retain the CPUE figures in the Sweave file.}
  \item{estM}{logical: if \code{TRUE} the routine will assume natural mortality \eqn{M}
    was estimated in the run, and the Sweave code uses the parameters \code{"M_1"} and 
    \code{"M_2"}; if \code{FALSE}, the function removes these parameters from the Sweave file.}
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
\note{
  Seave files can be found in the library directory: \cr
  \code{.../R/.../library/PBSawatea/snw/ymrrun-masterMCMC.Snw} \cr
  \code{.../R/.../library/PBSawatea/snw/ymrrun-masterMCMC20yrProjs.Snw} 
}
\seealso{
  \code{\link{runSweaveMCMC}}, \code{\link{runSweave}}, \code{\link{runMPD}}
}
\keyword{character}
\keyword{manip}


\name{runSweaveMCMC}
\alias{runSweaveMCMC}
\title{
  Run Customised Sweave Files for Awatea MCMC Runs
}
\description{
  Create and run customised Sweave files for \code{Awatea} MCMC runs.
}
\usage{
runSweaveMCMC( wd=getwd(), cpue=FALSE, estM=TRUE, 
    strSpp="YMR", filename="input25-ymr.txt", 
    runNo=25, rwtNo=0, running.awatea=0 )
}
\arguments{
  \item{wd}{working directory in which Awatea input and master Sweave files occur.}
  \item{cpue}{logical: if \code{TRUE} the routine will leave calls to 
    CPUE figures in the Sweave file, otherwise, they are removed.}
  \item{estM}{logical: if \code{TRUE} the routine will assume natural mortality \eqn{M}
    was estimated in the run, and the Sweave code uses the parameters \code{"M_1"} and 
    \code{"M_2"}; if \code{FALSE}, the function removes these parameters from the Sweave file.}
  \item{strSpp}{three-letter code that identifies the species.}
  \item{filename}{name of Awatea input file.}
  \item{runNo}{the run number that identifies a unique combination of input values.}
  \item{rwtNo}{the reweight number of the MPD run to build and collate figures.}
  \item{running.awatea}{numeric: if \code{0} load a previous \code{.rep} file, if \code{1} re-run \code{Awatea}}
}
\details{
  The values specified by the arguments (or derived variables) are directly 
  substituted into the Sweave file wherever similarly named variables 
  preceded by the \code{@} symbol occur.
  
  The \code{cpue} switch signals the removal of pieces of Sweave code if \code{cpue=TRUE}.
  
  The \code{estM} switch signals the removal of \code{"M_1"} and \code{"M_2"} from 
  the Sweave code if \code{estM=FALSE}.
}
\value{
  A customised Sweave file for \code{runNo} and \code{rwtNo} is created in a subdirectory \cr
  \code{./\emph{strSpp}run/MCMC.\emph{runNo}.\emph{rwtNo}} \cr
  from where it is run.
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
  \code{\link{runMCMC}}, \code{\link{runSweave}}, \code{\link{runADMB}}
}
\keyword{character}
\keyword{manip}


\name{runSweave}
\alias{runSweave}
\title{
  Run Customised Sweave Files for Awatea MPD Runs
}
\description{
  Create and run customised Sweave files for \code{Awatea} MPD runs.
}
\usage{
runSweave(wd=getwd(), cpue=FALSE, strSpp="XYZ",
    filename="spp-area-00.txt", runNo=1, rwtNo=0, running.awatea=0,
    Nsurvey=5, Snames=paste("Ser",1:Nsurvey,sep=""), delim="-")
}
\arguments{
  \item{wd}{working directory in which Awatea input and master Sweave files occur.}
  \item{cpue}{logical: if \code{TRUE} the routine will leave calls to 
    CPUE figures in the Sweave file, otherwise, they are removed.}
  \item{strSpp}{three-letter code that identifies the species.}
  \item{filename}{name of Awatea input file.}
  \item{runNo}{the run number that identifies a unique combination of input values.}
  \item{rwtNo}{the reweight number of the MPD run to build and collate figures.}
  \item{running.awatea}{numeric: if \code{0} load a previous \code{.rep} file, if \code{1} re-run \code{Awatea}}
  \item{Nsurvey}{number of surveys specified in the input file; controls image creation and placeholders.}
  \item{Snames}{survey names for figure labelling and image creation.}
  \item{delim}{character that delimits the components of the \code{filename} prefix.}
}
\details{
  The values specified by the arguments (or derived variables) are directly 
  substituted into the Sweave file wherever similarly named variables 
  preceded by the \code{@} symbol occur.
  
  The \code{cpue} switch signals the removal of pieces of Sweave code if \code{cpue=TRUE}.
  
  The \code{Nsurvey} argument essentially copies one Sweave line into \code{Nsurvey}
  lines for a set of lines identified by unique Sweave snippets.
}
\value{
  A customised Sweave file for \code{runNo} and \code{rwtNo} is created in a subdirectory \cr
  (e.g., \code{./YMRrun/MPD.29.01/}) from where it is run.
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
  \code{\link[PBSawatea]{runMPD}}, \code{\link[PBSawatea]{runSweaveMCMC}}, \code{\link[PBSawatea]{runADMB}}
}
\keyword{character}
\keyword{manip}


\name{runSweave}
\alias{runSweave}
\title{
  Run Customised Sweave Files for Awatea MPD Runs
}
\description{
  Create and run customised Sweave files for \code{Awatea} MPD runs.
}
\usage{
runSweave( wd = getwd(), strSpp = "XYZ",
   filename = "spp-area-00.txt", # Name of Awatea .txt file in 'run.dir' to run
   runNo   = 1,
   rwtNo   = 0,
   running.awatea = 0,          # 0 if just loading previous '.rep'; 1 if rerunning Awatea
   Nsex    = 2,                 # if 1 then Unisex, if 2 then Males & Females
   Ncpue   = 0,
   Nsurvey = 3,
   Snames  = paste("Ser",1:Nsurvey,sep=""),  # survey names (w/out spaces)
   SApos   = rep(TRUE,Nsurvey), # surveys with age composition data
   delim   = "-",
   debug   = FALSE,
   locode  = FALSE,             # source this function as local code (for development)
   awateaPath = "C:/Users/haighr/Files/Projects/ADMB/Coleraine",
   codePath   = "C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
   sexlab  = c("Females","Males")
) 
}
\arguments{
  \item{wd}{\emph{character}: working directory in which Awatea input files occur.}
  \item{strSpp}{\emph{character}: three-letter code that identifies the species.}
  \item{filename}{\emph{character}: name of Awatea input file.}
  \item{runNo}{\emph{numeric}: the run number that identifies a unique combination of input values.}
  \item{rwtNo}{\emph{numeric}: the reweight number of the MPD run to build and collate figures.}
  \item{running.awatea}{\emph{numeric}: if \code{0} load a previous \code{.rep} file, if \code{1} re-run \code{Awatea}}
  \item{Nsex}{\emph{numeric}: number of sexes; if \eqn{>1} then the routine will 
    activate Sweave lines with the prefix \code{@rmsex}, otherwise the lines are removed.}
  \item{Ncpue}{\emph{numeric}: number of CPUE series; if \eqn{>=1} then the routine will 
    activate Sweave lines with the prefix \code{@rmcpue}, otherwise the lines are removed.}
  \item{Nsurvey}{\emph{numeric}: number of survey series.}
  \item{Snames}{\emph{character}: survey names for figure labelling and image creation.}
  \item{SApos}{\emph{logical} vector of length \code{Nsurvey} denoting which surveys have ageing data.}
  \item{delim}{\emph{character} that delimits the components of the \code{filename} prefix.}
  \item{debug}{\emph{logical}: if \code{TRUE}, function only builds the target \code{.Snw} file but does not execute it.}
  \item{locode}{\emph{logical}: if \code{TRUE}, source local code located at \code{codePath}.}
  \item{awateaPath}{\emph{character}: system path to the \code{Awatea.exe} executable (using path delimiter `/`).}
  \item{codePath}{\emph{character}: system path to \pkg{PBSawatea} source code (using path delimiter `/`).}
  \item{sexlab}{\emph{character}: labels for the sexes in figures.}
}
\details{
  The values specified by the arguments (or derived variables) are directly 
  substituted into the Sweave file wherever similarly named variables 
  preceded by the \code{@} symbol occur.
  
  Additionally, if Sweave lines are preceded by the following control variables (prefixed \sQuote{\code{@rm}}),
  the lines are either activated or removed depending on function argument inputs.
  \tabular{cccc}{
    \bold{Sweave control} \tab \bold{function argument} \tab \bold{activate} \tab \bold{remove} \cr
    \code{@rmsex}    \tab \code{Nsex}   \tab \code{>1}   \tab \code{1} \cr
    \code{@rmcpue}   \tab \code{Ncpue}  \tab \code{>0}   \tab \code{0} 
  }
}
\value{
  A customised Sweave file for \code{runNo} and \code{rwtNo} is created in a subdirectory \cr
  (e.g., \code{./PSARC13/SGR/Data/Awatea/CST/SGRrun11/MPD.11.03}) from where it is run.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Example Sweave files (\code{run-master.Snw}) can be found in the library directory: \cr
  \code{.../R/.../library/PBSawatea/snw/}
}
\seealso{
  \code{\link[PBSawatea]{runMPD}}, \code{\link[PBSawatea]{runSweaveMCMC}}, \code{\link[PBSawatea]{runADMB}}
}
\keyword{character}
\keyword{manip}


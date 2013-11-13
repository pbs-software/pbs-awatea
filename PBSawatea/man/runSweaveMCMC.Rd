\name{runSweaveMCMC}
\alias{runSweaveMCMC}
\title{
  Run Customised Sweave Files for Awatea MCMC Runs
}
\description{
  Create and run customised Sweave files for \code{Awatea} MCMC runs.
}
\usage{
runSweaveMCMC (wd=getwd(), strSpp="XYZ",
   filename="spp-area-00.txt",  # name of Awatea .txt file in 'run.dir' to run
   runNo   = 1,
   rwtNo   = 0,
   running.awatea=0,            # if 0 then load prev. '.rep'; if 1 then rerun Awatea
   Nsex    = 2,                 # if 1 then Unisex, if 2 then Females & Males
   Ncpue   = 0,
   Nsurvey = 3,
   SApos   = rep(TRUE,Nsurvey), # surveys with age composition data
   mcsub   = 1:1000,
   delim   = "-",
   locode  = FALSE,             # source function as local code (for development)
   awateaPath = "C:/Users/haighr/Files/Projects/ADMB/Coleraine",
   codePath   = "C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
   histRP  = FALSE,             # historical reference points
   wpaper  = FALSE,             # working paper
   resdoc  = FALSE,             # research document
   redo.Graphs = TRUE           # recreate all the figures (.eps, .wmf, .png)
)
}
\arguments{
  \item{wd}{\emph{character}: working directory in which Awatea input files occur.}
  \item{strSpp}{\emph{character}: three-letter code that identifies the species.}
  \item{filename}{\emph{character}: name of Awatea input file.}
  \item{runNo}{\emph{numeric}: the run number that identifies a unique combination of input values.}
  \item{rwtNo}{\emph{numeric}: the reweight number of the MPD run to build and collate figures.}
  \item{running.awatea}{\emph{numeric}:  if \code{0} load a previous \code{.rep} file, if \code{1} re-run \code{Awatea}}
  \item{Nsex}{\emph{numeric}: number of sexes; if \eqn{>1} then the routine will 
    activate Sweave lines with the prefix \code{@rmsex}, otherwise the lines are removed.}
  \item{Ncpue}{\emph{numeric}: number of CPUE series; if \eqn{>=1} then the routine will 
    activate Sweave lines with the prefix \code{@rmcpue}, otherwise the lines are removed.}
  \item{Nsurvey}{\emph{numeric}: number of survey series.}
  \item{SApos}{\emph{logical} vector of length \code{Nsurvey} denoting which surveys have ageing data.}
  \item{mcsub}{\emph{numeric} vector of subsamples to select from the MCMC sample chain.}
  \item{delim}{\emph{character} that delimits the components of the \code{filename} prefix.}
  \item{locode}{\emph{logical}: if \code{TRUE}, source local code located at \code{codePath}.}
  \item{awateaPath}{\emph{character}: system path to the \code{Awatea.exe} executable (using path delimiter `/`).}
  \item{codePath}{\emph{character}: system path to \pkg{PBSawatea} source code (using path delimiter `/`).}
  \item{histRP}{\emph{logical}: if \code{TRUE}, the routine will activate Sweave lines with the 
    prefix \code{@rmhrp} (denoting historical reference points), otherwise the lines are removed.}
  \item{wpaper, resdoc}{\emph{logical}: if \code{TRUE}, the routine will remove Sweave lines with the 
    prefix \code{@rmresdoc} (denoting either working paper or research document), otherwise the lines are activated.
    Note that both \code{wpaper} and \code{resdoc} change the output file name.}
  \item{redo.Graphs}{\emph{logical}: if \code{TRUE}, recreate all the figures (\code{.eps}, \code{.wmf}, \code{.png}).}
}
\details{
  The values specified by the some of the arguments (and derived variables) are directly 
  substituted into the Sweave file wherever similarly named variables 
  preceded by the \code{@} symbol occur.
  
  Additionally, if Sweave lines are preceded by the following control variables (prefixed \sQuote{\code{@rm}}),
  the lines are either activated or removed depending on function argument inputs.
  \tabular{cccc}{
    \bold{Sweave control} \tab \bold{function argument} \tab \bold{activate} \tab \bold{remove} \cr
    \code{@rmsex}    \tab \code{Nsex}   \tab \code{>1}   \tab \code{1} \cr
    \code{@rmcpue}   \tab \code{Ncpue}  \tab \code{>0}   \tab \code{0} \cr
    \code{@rmhrp}    \tab \code{histRP} \tab \code{TRUE} \tab \code{FALSE} \cr
    \code{@rmresdoc} \tab \code{wpaper}, \code{resdoc} \tab \code{FALSE} \tab \code{TRUE} \cr
    \code{@rmROL}    \tab \code{strSpp} \tab \code{!="ROL"} \tab \code{=="ROL"}
  }

  The function \code{runSweaveMCMC} will incorporate various supplementary tex files (history and results) if they are available in the working directory.
  The results files must be named to match the run number (\code{runNo}) if they are to be incorporated into the MCMC tex document,
  and are substituted in the \code{run-masterMCMC.tex} build wherever the placeholder fields
  \code{@resultsMPD}, \code{@resultsMPDfigs}, \code{@resultsMPDtabs}, and \code{@resultsMCMC} occur.
  Only the run number is used when naming the supplementary tex files (see table below)
  as the species and area already determine the Awatea project structure 
  and the reweight number is likely chosen by the user to represent the run results.
  As an example, MCMC results for a coastwide (CST) Silvergray Rockfish (SGR) stock assessment might include 
  any of the following files depending on whether the output is for an informal MCMC summary, a working paper (WP), or a 
  research document (RD):
  \tabular{lll}{
    \bold{LaTeX file (\code{.tex})} \tab \bold{Placeholder} \tab \bold{Contents} \cr
    \code{runHistory}           \tab ---                    \tab {Summary of all model runs for CST SGR.} \cr
    \code{resultsMPD-run16}     \tab \code{@resultsMPD}     \tab {Summary of MPD results for Run 16 (WP/RD).}  \cr
    \code{resultsMPDfigs-run16} \tab \code{@resultsMPDfigs} \tab {Selected MPD figures created by \code{run-Master.Snw}.} \cr
    \code{resultsMPDtabs-run16} \tab \code{@resultsMPDtabs} \tab {Selected MPD tables created by \code{run-Master.Snw}.} \cr
    \code{resultsMCMC-run16}    \tab \code{@resultsMCMC}    \tab {Summary of MCMC results for Run 16 (WP/RD).}
  }
}
\value{
  A customised Sweave file for \code{runNo} and \code{rwtNo} is created in a subdirectory \cr
  (e.g., \code{./PSARC13/SGR/Data/Awatea/CST/SGRrun11/MCMC.11.03}) from where it is run.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Example Sweave files (\code{run-masterMCMC.Snw}) can be found in the library directory: \cr
  \code{.../R/.../library/PBSawatea/snw/}
}
\seealso{
  \code{\link[PBSawatea]{runMCMC}}, \code{\link[PBSawatea]{runSweave}}, \code{\link[PBSawatea]{runADMB}}
}
\keyword{character}
\keyword{manip}


\name{runSweaveMCMC}
\alias{runSweaveMCMC}
\alias{runMCMC}
\title{
  Run Customised Sweave Files for Awatea MCMC Runs
}
\description{
  Create and run customised Sweave files for \code{Awatea} MCMC runs.
}
\usage{
runSweaveMCMC (wd=getwd(), strSpp="XYZ",
   filename="spp-area-00.txt",
   runNo   = 1,
   rwtNo   = 0,
   running.awatea=0,
   Nsex    = 2,
   Ncpue   = 0,
   Nsurvey = 3,
   Ngear   = 1,
   Snames  = paste0("Ser",1:Nsurvey),
   SApos   = rep(TRUE,Nsurvey),
   Cnames  = paste0("Gear",1:Ngear),
   CApos   = rep(TRUE,Ngear),
   mcsub   = 1:1000,
   delim   = "-",
   locode  = FALSE,
   codePath=
     "C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
   histRP  = FALSE,
   wpaper  = FALSE,
   resdoc  = FALSE,
   redo.Graphs = TRUE,
   skip.last.year = TRUE,
   ptype   = "png",
   domeS   = FALSE,
   lang    = c("e","f")
)
runMCMC(prefix=c("spp","area"), runs=1, rwts=0, ...)
}
\arguments{
  \item{wd}{\code{character} -- working directory in which Awatea input files occur.}
  \item{strSpp}{\code{character} -- three-letter code that identifies the species.}
  \item{filename}{\code{character} -- name of Awatea input file.}
  \item{runNo}{\code{numeric} -- the run number that identifies a unique combination of input values.}
  \item{rwtNo}{\code{numeric} -- the reweight number of the MPD run to build and collate figures.}
  \item{running.awatea}{\code{numeric} -- if \code{0} load a previous \code{.rep} file, if \code{1} re-run \code{Awatea}}
  \item{Nsex}{\code{numeric} -- number of sexes; if \eqn{>1} then the routine will 
    activate Sweave lines with the prefix \code{@rmsex}, otherwise the lines are removed.}
  \item{Ncpue}{\code{numeric} -- number of CPUE series; if \eqn{>=1} then the routine will 
    activate Sweave lines with the prefix \code{@rmcpue}, otherwise the lines are removed.}
  \item{Nsurvey}{\code{numeric} -- number of survey series.}
  \item{Ngear}{\code{numeric} -- number of commercial gear series.}
  \item{Snames}{\code{character} -- survey names for figure labelling and image creation.}
  \item{SApos}{\code{logical} -- vector of length \code{Nsurvey} denoting which surveys have ageing data.}
  \item{Cnames}{\code{character} -- gear names for figure labelling and image creation.}
  \item{CApos}{\code{logical} -- vector of length \code{Ngear} denoting which commercial gears have ageing data.}
  \item{mcsub}{\code{numeric} -- vector of subsamples to select from the MCMC sample chain.}
  \item{delim}{\code{character} -- that delimits the components of the \code{filename} prefix.}
  \item{locode}{\code{logical} -- if \code{TRUE}, source local code located at \code{codePath}.}
  \item{codePath}{\code{character} -- system path to \pkg{PBSawatea} source code (using path delimiter `/`).}
  \item{histRP}{\code{logical} -- if \code{TRUE}, the routine will activate Sweave lines with the 
    prefix \code{@rmhrp} (denoting historical reference points), otherwise the lines are removed.}
  \item{wpaper, resdoc}{\code{logical} -- if \code{TRUE}, the routine will remove Sweave lines with the 
    prefix \code{@rmresdoc} (denoting either working paper or research document), otherwise the lines are activated.
    Note that both \code{wpaper} and \code{resdoc} change the output file name.}
  \item{redo.Graphs}{\code{logical} -- if \code{TRUE}, recreate all the figures (\code{.eps}, \code{.wmf}, \code{.png}).}
  \item{skip.last.year}{\code{logical} -- if \code{TRUE}, remove the last year of projections
    (set to FALSE for POP 5ABC in 2010)}
  \item{ptype}{\code{character} -- plot type for the figures: \code{"png"} = portable network graphics, 
    \code{"eps"} = encapsulated postscript.}
  \item{domeS}{\code{logical} -- if \code{TRUE}, estimate the right-hand side of the selectivity curve,
    i.e., assume dome-shaped selectivity.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).\cr\cr------ \code{runMCMC} ------}


  \item{prefix}{\code{character} -- vector of \code{filename} prefix components delimited by \code{delim}.}
  \item{runs}{\code{numeric} -- vector of run number(s) for \code{runNo} in \code{runSweaveMCMC}.}
  \item{rwts}{\code{numeric} -- vector of reweight number(s) for \code{rwtNo} in \code{runSweaveMCMC}.}
  \item{...}{arguments passed on to the function \code{runSweaveMCMC}.}
}
\details{
  The values specified by the some of the arguments (and derived variables) are directly 
  substituted into the Sweave file \code{run-masterMCMC.Snw} wherever similarly named variables 
  preceded by the \code{@} symbol occur.
  Some primary substitutions when building the temporary Sweave file (\code{tfile}):
  \preformatted{
   tfile = gsub("@cwd",wd,tfile)
   tfile = gsub("@model.name",model.name,tfile)
   tfile = gsub("@run.dir",run.dir,tfile)
   tfile = gsub("@fig.dir",mc.dir,tfile)
   tfile = gsub("@mpd.dir",mpd.dir,tfile)
   tfile = gsub("@msy.dir",msy.dir,tfile)
   tfile = gsub("@prj.dir",prj.dir,tfile)
   tfile = gsub("@running.awatea",running.awatea,tfile)
   tfile = gsub("@redo.Graphs",redo.Graphs,tfile)
   tfile = gsub("@skip.last.year",skip.last.year,tfile)
   tfile = gsub("@mcsub",deparse(mcsub),tfile)
   tfile = gsub("@nsex",Nsex,tfile)
   tfile = gsub("@ngear",Ngear,tfile)
   tfile = gsub("@sppcode",strSpp,tfile)
   tfile = gsub("@ptype",ptype,tfile)
   tfile = gsub("@lang",deparse(lang),tfile)
   }
  Additionally, if Sweave lines are preceded by the following control variables (prefixed \sQuote{\code{@rm}}),
  the lines are either activated or removed depending on function argument inputs.
  \tabular{cccc}{
    \bold{Sweave control} \tab \bold{function argument} \tab \bold{activate when} \tab \bold{remove when} \cr
    \code{@rmsex}    \tab \code{Nsex}          \tab \code{>1}    \tab \code{1} \cr
    \code{@rmcpue}   \tab \code{Ncpue}         \tab \code{>0}    \tab \code{0} \cr
    \code{@rmCA}     \tab \code{sum(CApos)}    \tab \code{>0}    \tab \code{0} \cr
    \code{@rmSA}     \tab \code{sum(SApos)}    \tab \code{>0}    \tab \code{0} \cr
    \code{@rmhrp}    \tab \code{histRP}        \tab \code{TRUE}  \tab \code{FALSE} \cr
    \code{@rmresdoc} \tab \code{wpaper||resdoc}\tab \code{FALSE} \tab \code{TRUE} \cr
    \code{@rmROL}    \tab \code{strSpp=="ROL"} \tab \code{FALSE} \tab \code{TRUE}\cr
    \code{@rmdome}   \tab \code{domeS}         \tab \code{TRUE}  \tab \code{FALSE}
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
    \code{resultsMPDfigs-run16} \tab \code{@resultsMPDfigs} \tab {Selected MPD figures created by \code{run-master.Snw}.} \cr
    \code{resultsMPDtabs-run16} \tab \code{@resultsMPDtabs} \tab {Selected MPD tables created by \code{run-master.Snw}.} \cr
    \code{resultsMCMC-run16}    \tab \code{@resultsMCMC}    \tab {Summary of MCMC results for Run 16 (WP/RD).}
  }

  The wrapper function \code{runMCMC} essentially loops through \code{runSweaveMCMC} 
  using \code{(i in runs)} and \code{(j in rwts)}.
}
\value{
  A customised Sweave file for \code{runNo} and \code{rwtNo} is created in a subdirectory \cr
  (e.g., \code{./PSARC13/SGR/Data/Awatea/CST/SGRrun11/MCMC.11.03}) from where it is run.

  \code{runMCMC} produces multiple Sweave documents and PDF files.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-25}
}
\note{
  Example Sweave files (\code{run-masterMCMC.Snw}) can be found in the library directory: \cr
  \code{paste0(system.file(package="PBSawatea"),"/snw")}
}
\seealso{
  \code{\link[PBSawatea]{runSweave}}, \code{\link[PBSawatea]{runADMB}}
}
\keyword{character}
\keyword{manip}


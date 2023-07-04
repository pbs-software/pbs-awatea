\name{runSweave}
\alias{runSweave}
\alias{runSweaveMCMC}
\alias{runMPD}
\alias{runMCMC}
\title{
  Run Customised Sweave Files for Awatea MPD and MCMC Runs
}
\description{
  Create and run customised Sweave files for \code{Awatea} MPD and MCMC runs.
}
\usage{
runSweave(wd=getwd(), strSpp="XYZ", filename="spp-area-00.txt",
   running.awatea=0, runNo=1, rwtNo=0, Nsex=2, Ncpue=0, Nsurvey=3,
   Ngear=1, NCAset=1, maxcol=5, sexlab=c("Females","Males"), 
   Snames=paste0("Ser", 1:Nsurvey), SApos=rep(TRUE, Nsurvey), 
   Cnames=paste0("Gear", 1:Ngear),  CApos=rep(TRUE, Ngear),
   Unames=paste0("CPUE", 1:Ncpue), delim="-", debug=FALSE, 
   dome=FALSE, locode=FALSE, 
   codePath="C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
   resdoc=FALSE, redo.Graphs=TRUE, ptype="png", lang=c("e", "f")) 

runSweaveMCMC(wd=getwd(), strSpp="XYZ", filename="spp-area-00.txt",
   running.awatea=0, runNo=1, rwtNo=0, Nsex=2, Ncpue=0, Nsurvey=3,
   Ngear=1, Snames=paste0("Ser", 1:Nsurvey), SApos=rep(TRUE, Nsurvey), 
   Cnames=paste0("Gear", 1:Ngear), CApos=rep(TRUE, Ngear),
   Unames=paste0("CPUE", 1:Ncpue), delim="-", dome=FALSE, locode=FALSE,
   codePath="C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
   resdoc=FALSE, redo.Graphs=TRUE, ptype="png", lang=c("e", "f"),
   mcsub=1:1000, histRP=FALSE, skip.last.year=TRUE)

runMPD(prefix=c("spp","area"), runs=1, rwts=0, ...)
runMCMC(prefix=c("spp","area"), runs=1, rwts=0, ...)
}
\arguments{
  \item{wd}{\code{character} -- working directory in which Awatea input files occur.}
  \item{strSpp}{\code{character} -- three-letter code that identifies the species.}
  \item{filename}{\code{character} -- name of Awatea input file.}
  \item{running.awatea}{\code{numeric} -- if \code{0} load a previous \code{.rep} file, if \code{1} re-run \code{Awatea}}
  \item{runNo}{\code{numeric} -- the run number that identifies a unique combination of input values.}
  \item{rwtNo}{\code{numeric} -- the reweight number of the MPD run to build and collate figures.}
  \item{Nsex}{\code{numeric} -- number of sexes; if \eqn{>1} then the routine will 
    activate Sweave lines with the prefix \code{@rmsex}, otherwise the lines are removed.}
  \item{Ncpue}{\code{numeric} -- number of CPUE series; if \eqn{>=1} then the routine will 
    activate Sweave lines with the prefix \code{@rmcpue}, otherwise the lines are removed.}
  \item{Nsurvey}{\code{numeric} -- number of survey series.}
  \item{Ngear}{\code{numeric} -- number of commercial gear series.}
  \item{NCAset}{\code{numeric} -- number of commercial catch-age-age plot sets (separate figures).
    It is best to limit the number of panels (years) per page to 20; therefore 
    \code{NCAset=1} for \eqn{\le}{<=}20 CA years, \code{NCAset=2} for \eqn{\le}{<=}40 CA years, and so on.}
  \item{maxcol}{\code{numeric} -- maximum number of columns per page for CA plots.}
  \item{sexlab}{\code{character} -- labels for the sexes in figures.}
  \item{Snames}{\code{character} -- survey names for figure labelling and image creation.}
  \item{SApos}{\code{logical} -- vector of length \code{Nsurvey} denoting which surveys have ageing data.}
  \item{Cnames}{\code{character} -- commercial fishery names for figure labelling and image creation.}
  \item{CApos}{\code{logical} -- vector of length \code{Ngear} denoting which commercial gears have ageing data.}
  \item{Unames}{\code{character} -- CPUE names for figure labelling}
  \item{delim}{\code{character} -- that delimits the components of the \code{filename} prefix.}
  \item{debug}{\code{logical} -- if \code{TRUE}, function only builds the target \code{.Snw} file but does not execute it.}
  \item{dome}{\code{logical} -- if \code{TRUE}, the variance of the descending limb of the selectivity curve has been estimated}
  \item{locode}{\code{logical} -- if \code{TRUE}, source local code located at \code{codePath}.}
  \item{codePath}{\code{character} -- system path to \pkg{PBSawatea} source code (using path delimiter `/`).}
  \item{resdoc}{\code{logical} -- if \code{TRUE}, the routine will remove Sweave lines with the 
    prefix \code{@rmresdoc} (denoting either working paper or research document), otherwise the lines are activated.}
  \item{redo.Graphs}{\code{logical} -- if \code{TRUE}, re-create the figures.}
  \item{ptype}{\code{character} -- plot type for the figures: \code{"png"} = portable network graphics, 
    \code{"eps"} = encapsulated postscript.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
  \item{mcsub}{\code{numeric} -- indices to subset the saved samples from MCMC simulations to remove burn-ins}
  \item{histRP}{\code{logical} -- if \code{TRUE}, the routine will activate Sweave lines with the 
    prefix \code{@rmhrp} (denoting historical reference points), otherwise the lines are removed.}
  \item{skip.last.year}{\code{logical} -- if \code{TRUE}, remove last year of projections (set to \code{FALSE} for POP 5ABC in 2010)\cr\cr------ \code{runMPD|runMCMC} ------}

  \item{prefix}{\code{character} -- vector of \code{filename} prefix components delimited by \code{delim}.}
  \item{runs}{\code{numeric} -- vector of run number(s) for \code{runNo} in \code{runSweave|runSweaveMCMC}.}
  \item{rwts}{\code{numeric} -- vector of reweight number(s) for \code{rwtNo} in \code{runSweave|runSweaveMCMC}.}
  \item{...}{\code{dots} -- additional arguments passed to the function \code{runSweave|runSweaveMCMC}.}
}
\details{
  The values specified by the arguments (or derived variables) are directly 
  substituted into the Sweave file \code{run-master.Snw|run-masterMCMC.Snw} wherever similarly named variables 
  preceded by the \code{@} symbol occur. 
  Some primary substitutions when building the temporary Sweave file (\code{tfile}) for MPD results:
  \preformatted{
   tfile = gsub("@cwd",wd,tfile)
   tfile = gsub("@model.name",model.name,tfile)
   tfile = gsub("@run.dir",run.dir,tfile)
   tfile = gsub("@fig.dir",mpd.dir,tfile)
   tfile = gsub("@running.awatea",running.awatea,tfile)
   tfile = gsub("@redo.Graphs",redo.Graphs,tfile)
   tfile = gsub("@sexlab",deparse(sexlab),tfile)
   tfile = gsub("@sppcode",strSpp,tfile)
   tfile = gsub("@ptype",ptype,tfile)
   tfile = gsub("@lang",deparse(lang),tfile)
   }
  Additional substitutions when building the temporary Sweave file (\code{tfile}) for MCMC results:
  \preformatted{
   tfile = gsub("@mpd.dir",mpd.dir,tfile)
   tfile = gsub("@msy.dir",msy.dir,tfile)
   tfile = gsub("@prj.dir",prj.dir,tfile)
   tfile = gsub("@nsex",Nsex,tfile)
   tfile = gsub("@ngear",Ngear,tfile)
   tfile = gsub("@mcsub",deparse(mcsub),tfile)
   tfile = gsub("@skip.last.year",skip.last.year,tfile)
   }
  Furthermore, if Sweave lines are preceded by the following control variables (prefixed \sQuote{\code{@rm}}),
  the lines are either activated or removed depending on function argument inputs.
  \tabular{cccc}{
    \bold{Sweave control} \tab \bold{function argument} \tab \bold{activate when} \tab \bold{remove when}\cr
    \code{@rmsex}    \tab \code{Nsex}          \tab \code{>1}    \tab \code{1}\cr
    \code{@rmcpue}   \tab \code{Ncpue}         \tab \code{>0}    \tab \code{0}\cr
    \code{@rmCA}     \tab \code{sum(CApos)}    \tab \code{>0}    \tab \code{0}\cr
    \code{@rmSA}     \tab \code{sum(SApos)}    \tab \code{>0}    \tab \code{0}\cr
    \code{@rmCSA}    \tab \code{sum(CApos)+sum(SApos)} \tab \code{>0} \tab \code{0}\cr
    \code{@rmCA1}    \tab \code{NCAset}        \tab \code{1}     \tab \code{!=1} \cr
    \code{@rmCA2}    \tab \code{NCAset}        \tab \code{2}     \tab \code{!=2} \cr
    \code{@rmresdoc} \tab \code{resdoc}        \tab \code{FALSE} \tab \code{TRUE}\cr
    \code{@rmROL}    \tab \code{strSpp=="ROL"} \tab \code{FALSE} \tab \code{TRUE}\cr
    \code{@rmdome}   \tab \code{domeS}         \tab \code{TRUE}  \tab \code{FALSE}\cr
    \code{@rmhrp}    \tab \code{histRP}        \tab \code{TRUE}  \tab \code{FALSE}
  }
  The wrapper functions \code{runMPD|runMCMC} essentially loop through \code{runSweave|runSweaveMCMC},
  respectively, using \code{(i in runs)} and \code{(j in rwts)}.
}
\value{
  A customised Sweave file for \code{runNo} and \code{rwtNo} is created in a subdirectory \cr
  (e.g., \code{./PSARC13/SGR/Data/Awatea/CST/SGRrun11/MPD.11.03}) from where it is run.

  \code{runMPD|runMCMC} produce multiple Sweave documents and PDF files.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Regional Headquarters (RHQ), Vancouver BC\cr
  Last modified \code{Rd: 2023-06-29}
}
\note{
  Example Sweave files (\code{run-master.Snw|run-masterMCMC.Snw}) can be found in the library directory: \cr
  \code{system.file("snw",package="PBSawatea")}
}
\seealso{
\pkg{PBSawatea}:
  \code{\link[PBSawatea]{runADMB}}
}
\keyword{character}
\keyword{manip}


\name{runSweave}
\alias{runSweave}
\alias{runMPD}
\title{
  Run Customised Sweave Files for Awatea MPD Runs
}
\description{
  Create and run customised Sweave files for \code{Awatea} MPD runs.
}
\usage{
runSweave( wd = getwd(), strSpp="XYZ",
   filename = "spp-area-00.txt",
   runNo   = 1,
   rwtNo   = 0,
   running.awatea =0,
   Nsex    = 2,
   Ncpue   = 0,
   Nsurvey = 3,
   Ngear   = 1,
   NCAset  = 1,
   Snames  = paste0("Ser",1:Nsurvey),
   SApos   = rep(TRUE,Nsurvey),
   Cnames  = paste0("Gear",1:Ngear),
   CApos   = rep(TRUE,Ngear),
   delim   = "-",
   debug   = FALSE,
   locode  = FALSE,
   codePath=
     "C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
   sexlab  = c("Females","Males"),
   resdoc  = FALSE,
   redo.Graphs = TRUE,
   ptype   = "png",
   lang    = c("e","f")
) 
runMPD(prefix=c("spp","area"), runs=1, rwts=0, ...)
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
  \item{NCAset}{\code{numeric} -- number of commercial catch-age-age plot sets (separate figures).
    It is best to limit the number of panels (years) per page to 20; therefore 
    \code{NCAset=1} for \eqn{\le}{<=}20 CA years, \code{NCAset=2} for \eqn{\le}{<=}40 CA years, and so on.}
  \item{Snames}{\code{character} -- survey names for figure labelling and image creation.}
  \item{SApos}{\code{logical} -- vector of length \code{Nsurvey} denoting which surveys have ageing data.}
  \item{Cnames}{\code{character} -- gear names for figure labelling and image creation.}
  \item{CApos}{\code{logical} -- vector of length \code{Ngear} denoting which commercial gears have ageing data.}
  \item{delim}{\code{character} -- that delimits the components of the \code{filename} prefix.}
  \item{debug}{\code{logical} -- if \code{TRUE}, function only builds the target \code{.Snw} file but does not execute it.}
  \item{locode}{\code{logical} -- if \code{TRUE}, source local code located at \code{codePath}.}
  \item{codePath}{\code{character} -- system path to \pkg{PBSawatea} source code (using path delimiter `/`).}
  \item{sexlab}{\code{character} -- labels for the sexes in figures.}
  \item{resdoc}{\code{logical} -- if \code{TRUE}, the routine will remove Sweave lines with the 
    prefix \code{@rmresdoc} (denoting either working paper or research document), otherwise the lines are activated.}
  \item{redo.Graphs}{\code{logical} -- if \code{TRUE}, re-create the figures.}
  \item{ptype}{\code{character} -- plot type for the figures: \code{"png"} = portable network graphics, 
    \code{"eps"} = encapsulated postscript.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).\cr\cr------ \code{runMPD} ------}

  \item{prefix}{\code{character} -- vector of \code{filename} prefix components delimited by \code{delim}.}
  \item{runs}{\code{numeric} -- vector of run number(s) for \code{runNo} in \code{runSweave}.}
  \item{rwts}{\code{numeric} -- vector of reweight number(s) for \code{rwtNo} in \code{runSweave}.}
  \item{...}{arguments passed on to the function \code{runSweave}.}
}
\details{
  The values specified by the arguments (or derived variables) are directly 
  substituted into the Sweave file \code{run-master.Snw} wherever similarly named variables 
  preceded by the \code{@} symbol occur. 
  Some primary substitutions when building the temporary Sweave file (\code{tfile}):
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
  Additionally, if Sweave lines are preceded by the following control variables (prefixed \sQuote{\code{@rm}}),
  the lines are either activated or removed depending on function argument inputs.
  \tabular{cccc}{
    \bold{Sweave control} \tab \bold{function argument} \tab \bold{activate when} \tab \bold{remove when} \cr
    \code{@rmsex}    \tab \code{Nsex}          \tab \code{>1}    \tab \code{1} \cr
    \code{@rmcpue}   \tab \code{Ncpue}         \tab \code{>0}    \tab \code{0} \cr
    \code{@rmCA}     \tab \code{sum(CApos)}    \tab \code{>0}    \tab \code{0} \cr
    \code{@rmSA}     \tab \code{sum(SApos)}    \tab \code{>0}    \tab \code{0} \cr
    \code{@rmCSA}    \tab \code{sum(CApos)+sum(SApos)} \tab \code{>0} \tab \code{0} \cr
    \code{@rmCA1}    \tab \code{NCAset}        \tab \code{1}     \tab \code{!=1} \cr
    \code{@rmCA2}    \tab \code{NCAset}        \tab \code{2}     \tab \code{!=2} \cr
    \code{@rmresdoc} \tab \code{resdoc}        \tab \code{FALSE} \tab \code{TRUE} \cr
    \code{@rmROL}    \tab \code{strSpp=="ROL"} \tab \code{FALSE} \tab \code{TRUE}
  }
  The wrapper function \code{runMPD} essentially loops through \code{runSweave} 
  using \code{(i in runs)} and \code{(j in rwts)}.
}
\value{
  A customised Sweave file for \code{runNo} and \code{rwtNo} is created in a subdirectory \cr
  (e.g., \code{./PSARC13/SGR/Data/Awatea/CST/SGRrun11/MPD.11.03}) from where it is run.

  \code{runMPD} produces multiple Sweave documents and PDF files.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-25}
}
\note{
  Example Sweave files (\code{run-master.Snw}) can be found in the library directory: \cr
  \code{.../R/.../library/PBSawatea/snw/}
}
\seealso{
  \code{\link[PBSawatea]{runSweaveMCMC}}, \code{\link[PBSawatea]{runADMB}}
}
\keyword{character}
\keyword{manip}


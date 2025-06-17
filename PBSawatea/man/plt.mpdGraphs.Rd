\name{plt.mpdGraphs}
\alias{plt.mpdGraphs}
\alias{plt.mcmcGraphs}
\title{
  Plot a Set of Figures for MPD and MCMC Runs
}
\description{
  Plot a set of portable network graphics (\code{.png}) (and/or
  encapsulated postscript (\code{.eps})) figures for the MPD(mode of
  the posterior distribution) and MCMC (Monte Carlo Markoff Chain) results.
}
\usage{
plt.mpdGraphs(obj, save=FALSE, ssnames=paste("Ser",1:9,sep=""),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, ngear=1,
   pchGear=seq(21,20+ngear,1), ltyGear=seq(1,ngear,1), 
   colGear=rep(c("black","blue"),ngear)[1:ngear], lang=c("e","f"))

plt.mcmcGraphs(mcmcObj, projObj=NULL, mpdObj=NULL, save=FALSE, 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, ngear=1,
   ylim.recruitsMCMC=NULL, ylim.exploitMCMC=NULL,
   ylim.VBcatch=NULL, ylim.BVBnorm=NULL,
   xlim.snail=NULL, ylim.snail=NULL,
   plotPolicies=names(projObj$Y[1:6]),
   onePolicy=names(projObj$Y[2]), mpd=list(),
   SAR.width=7.5, SAR.height=4, trevObj=NULL, lang=c("e","f"))
}
\arguments{
  \item{obj}{\code{list} -- an Awatea MPD results object (\emph{e.g.}, \code{currentRes})}
  \item{save}{\code{logical} -- \bold{not used} \emph{deprecated}}
  \item{ssnames}{\code{character} -- survey series names for figure labelling and image creation}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{eps} = encapsulated postscript, \code{png} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{png} file.}
  \item{ngear}{\code{numeric} -- number of commercial gears}
  \item{pchGear, ltyGear, colGear}{\code{numeric|character} -- symbols, line types, and colours for plotting data by commercial gear.}
  \item{mcmcObj}{\code{list} -- an Awatea MCMC object (\emph{e.g.}, \code{currentMCMC})}
  \item{projObj}{\code{list} -- an Awatea projected biomass object (\emph{e.g.}, \code{currentProj})}
  \item{mpdObj}{\code{list} -- an Awatea MPD results object (\emph{e.g.}, \code{currentRes})}
  \item{ylim.recruitsMCMC}{\code{numeric} -- y-limits for the recruitsMCMC plot}
  \item{ylim.exploitMCMC}{\code{numeric} -- y-limits for the recruitsMCMC plot}
  \item{ylim.VBcatch}{\code{numeric} -- y-limits for the VBcatch plot}
  \item{ylim.BVBnorm}{\code{numeric} -- y-limits for the BVBnorm plot}
  \item{xlim.snail}{\code{numeric} -- x-limits for the snail-trail plot}
  \item{ylim.snail}{\code{numeric} -- y-limits for the snail-trail plot}
  \item{plotPolicies}{\code{numeric} -- six-policies projections to plot}
  \item{onePolicy}{\code{numeric} -- the one policy from the above six to use for some figures}
  \item{mpd}{\code{list} -- list of MPD values for parameters, biomass, and recruitment}
  \item{SAR.width}{\code{numeric} -- width (inches) of an output WMF file for SAR figures}
  \item{SAR.height}{\code{numeric} -- height (inches) of an output WMF file for SAR figures}
  \item{trevObj}{\code{data.frame} -- various MCMC values suggested by Trevor Davies (Dalhousie) for comparison:
    \code{c("h", "M_1", "M_2", "B0", "MSY", "Bmsy", "umsy")}. }
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
}
\details{
  Creates a whole heap of encapsulated postscript (\code{.eps}) and/or 
  portable network graphics (\code{.png}) files, depending on the plot 
  types specified in \code{ptypes}.
  Due to endless problems with postscript files (not to mention
  their potentially humungous sizes), \code{.png} files are favoured
  for stock assessment documents by the CSAP/CSAS office.

  As of July 2018, english plot are dumped into the current working
  directory (as per usual) while french plots, using the same file names,
  are dumped into a subdirectory called \sQuote{\code{french}}.
}
\value{
  \bold{plt.mpdGraphs} creates the following PNG (and/or EPS) files:
  \tabular{ll}{
    \code{exploit.png}          \tab annual exploitation rate\cr
    \code{recruits.png}         \tab annual recruitment at age 1\cr
    \code{selectivity.png}      \tab selectivity curves for commercial gear(s) and survey(s)\cr
    \code{ageComm.png}          \tab fits to annual commercial age composition (panel plots)\cr
    \code{ageSurv.png}          \tab fits to annual survey age composition (panel plots)\cr
    \code{survIndSer.png}       \tab four figures of survey indices (calls \code{\link{plotIndexNotLattice}})\cr
    \code{CPUEser.png....}      \tab CPUE indices with error bars (calls \code{\link{plotCPUE}})\cr
    \code{commAgeResids.png}    \tab standardised age residuals for commercial gear by sex\cr
    \code{survAgeResidsSer.png} \tab standardised age residuals for surveys by sex\cr
    \code{meanAge.png}          \tab mean age for catch and surveys\cr
    \code{stockRecruit.png}     \tab stock recruitment function
  }

  \bold{plt.mcmcGraphs} creates the following postscript files: \cr
  \tabular{ll}{
    \code{recruitsMCMC.png}    \tab boxplots of annual MCMC recruitment\cr
    \code{exploitMCMC.png}     \tab boxplots of annual MCMC exploitation rate\cr
    \code{pdfParameters.png}   \tab density panel plots of MCMC parameter priors and posteriors\cr
    \code{pdfBiomass.png}      \tab density panel plots of annual female spawning biomass\cr
    \code{pdfRecruitment.png}  \tab density panel plots of annual recuitment\cr
    \code{selectivityMCMC.png} \tab \bold{not currently implemented}\cr
    \code{traceBiomass.png}    \tab panel plots of annual spawning biomass with running quantiles\cr
    \code{traceRecruits.png}   \tab panel plots of annual recruitmwent traces with running quantiles\cr
    \code{traceParams.png}     \tab panel plots of parameter traces with running quantiles\cr
    \code{splitChain.png}      \tab panel plots of cumulative parameter estimate chains\cr
    \code{VBcatch.png}         \tab boxplots of annual vulnerable biomass and barplots of catch\cr
    \code{BVBnorm.png}         \tab spawning and vulnerable biomass relative to their virgin levels\cr
    \code{Bproj.png}           \tab boxplots of spawning biomass -- MCMCs and projections\cr
    \code{Rproj.png}           \tab boxplots of recruitment -- MCMCs and projections\cr
    \code{RprojOnePolicy.png}  \tab boxplots of recruitment for one catch policy\cr
    \code{snail.png}           \tab time series of \eqn{u_{t-1}/u_\mathrm{MSY}}{u[t-1]/u[MSY]} \emph{vs.} \eqn{B_t/B_\mathrm{MSY}}{B[t]/B[MSY]}\cr
    \code{pairs[1,...,n].png}  \tab pairs plot of parameter MCMC samples (6 pairs/page)\cr
    \code{pairsMSY.png}        \tab pairs plot of MCMC samples for MSY-related parameters
  }
}
\author{
First created by:\cr
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew M. Edwards}, PBS, DFO, Nanaimo BC

Modified by:\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-17}
}
\note{
  The function is currently set up to send figures to portable network graphics
  (png) and/or encapsulated postscript (eps) files only. If a user wants other file types
  (e.g., windows metafiles), then additional lines of code need to be added to invoke
  the appropriate device(s).
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{plt.mcmcGraphs}},
  \code{\link[PBSawatea]{panelChains}},
  \code{\link[PBSawatea]{plotIndexNotLattice}},
  \code{\link[PBSawatea]{plotCPUE}},
  \code{\link[PBSawatea]{plotSnail}}
}
\keyword{hplot}
\keyword{connection}% __ONLY ONE__ keyword per line
\concept{PBSscape.r}


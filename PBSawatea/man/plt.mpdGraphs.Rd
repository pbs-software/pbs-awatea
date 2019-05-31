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
  \bold{plt.mpdGraphs} creates the following PNG (and/or EPS) files: \cr
  \code{exploit.png............}{annual exploitation rate,} \cr
  \code{recruits.png...........}{annual recruitment at age 1,} \cr
  \code{selectivity.png........}{selectivity curves for commercial gear(s) and survey(s),} \cr
  \code{ageComm.png............}{fits to annual commercial age composition (panel plots),} \cr
  \code{ageSurv.png............}{fits to annual survey age composition (panel plots),} \cr
  \code{survIndSer.png.........}{four figures of survey indices (calls \code{\link{plotIndexNotLattice}}),} \cr
  \code{CPUEser.png............}{CPUE indices with error bars (calls \code{\link{plotCPUE}}),} \cr
  \code{commAgeResids.png......}{standardised age residuals for commercial gear by sex,} \cr
  \code{survAgeResidsSer.png...}{standardised age residuals for surveys by sex,} \cr
  \code{meanAge.png............}{mean age for catch and surveys,} \cr
  \code{stockRecruit.png.......}{stock recruitment function.} \cr

  \bold{plt.mcmcGraphs} creates the following postscript files: \cr
  \code{recruitsMCMC.png.......}{boxplots of annual MCMC recruitment,} \cr
  \code{exploitMCMC.png........}{boxplots of annual MCMC exploitation rate,} \cr
  \code{pdfParameters.png......}{density panel plots of MCMC parameter priors and posteriors,} \cr
  \code{pdfBiomass.png.........}{density panel plots of annual female spawning biomass,} \cr
  \code{pdfRecruitment.png.....}{density panel plots of annual recuitment,} \cr
  \code{selectivityMCMC.png....}{\bold{not currently implemented},} \cr
  \code{traceBiomass.png.......}{panel plots of annual spawning biomass with running quantiles,} \cr
  \code{traceRecruits.png......}{panel plots of annual recruitmwent traces with running quantiles,} \cr
  \code{traceParams.png........}{panel plots of parameter traces with running quantiles,} \cr
  \code{splitChain.png.........}{panel plots of cumulative parameter estimate chains,} \cr
  \code{VBcatch.png............}{boxplots of annual vulnerable biomass and barplots of catch,} \cr
  \code{BVBnorm.png............}{spawning and vulnerable biomass relative to their virgin levels,} \cr
  \code{Bproj.png..............}{boxplots of spawning biomass -- MCMCs and projections,} \cr
  \code{Rproj.png..............}{boxplots of recruitment -- MCMCs and projections,} \cr
  \code{RprojOnePolicy.png.....}{boxplots of recruitment for one catch policy,} \cr
  \code{snail.png..............}{time series of \eqn{u_{t-1}/u_\mathrm{MSY}}{u[t-1]/u[MSY]} \emph{vs.} \eqn{B_t/B_\mathrm{MSY}}{B[t]/B[MSY]},} \cr
  \code{pairs[1,...,n].png.....}{pairs plot of parameter MCMC samples (6 pairs/page).} \cr
  \code{pairsMSY.png...........}{pairs plot of MCMC samples for MSY-related parameters.}
}
\author{
First created by:\cr
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew M. Edwards}, PBS, DFO, Nanaimo BC

Modified by:\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2019-05-31}
}
\note{
  The function is currently set up to send figures to portable network graphics
  (png) and/or encapsulated postscript (eps) files only. If a user wants other file types
  (e.g., windows metafiles), then additional lines of code need to be added to invoke
  the appropriate device(s).
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link{plt.mcmcGraphs}},
  \code{\link{panelChains}},
  \code{\link{plotIndexNotLattice}},
  \code{\link{plotCPUE}},
  \code{\link{plotSnail}}
}
\keyword{hplot}
\keyword{connection}% __ONLY ONE__ keyword per line


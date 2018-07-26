\name{plt.mpdGraphs}
\alias{plt.mpdGraphs}
\alias{plt.mcmcGraphs}
\title{
  Plot a Set of Figures for MPD and MCMC Runs
}
\description{
  Plot a set of encapsulated postscript (\code{.eps}) and portable network graphics (\code{.png}) figures  for the MPD (mode of the posterior distribution)
  and MCMC (Monte Carlo Markoff Chain) results.
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
  \bold{plt.mpdGraphs} creates the following EPS (and/or PNG) files: \cr
  \code{exploit.eps............}{annual exploitation rate,} \cr
  \code{recruits.eps...........}{annual recruitment at age 1,} \cr
  \code{selectivity.eps........}{selectivity curves for commercial gear(s) and survey(s),} \cr
  \code{ageComm.eps............}{fits to annual commercial age composition (panel plots),} \cr
  \code{ageSurv.eps............}{fits to annual survey age composition (panel plots),} \cr
  \code{survIndSer.eps.........}{four figures of survey indices (calls \code{\link{plotIndexNotLattice}}),} \cr
  \code{CPUEser.eps............}{CPUE indices with error bars (calls \code{\link{plotCPUE}}),} \cr
  \code{commAgeResids.eps......}{standardised age residuals for commercial gear by sex,} \cr
  \code{survAgeResidsSer.eps...}{standardised age residuals for surveys by sex,} \cr
  \code{meanAge.eps............}{mean age for catch and surveys,} \cr
  \code{stockRecruit.eps.......}{stock recruitment function.} \cr

  \bold{plt.mcmcGraphs} creates the following postscript files: \cr
  \code{recruitsMCMC.eps.......}{boxplots of annual MCMC recruitment,} \cr
  \code{exploitMCMC.eps........}{boxplots of annual MCMC exploitation rate,} \cr
  \code{pdfParameters.eps......}{density panel plots of MCMC parameter priors and posteriors,} \cr
  \code{pdfBiomass.eps.........}{density panel plots of annual female spawning biomass,} \cr
  \code{pdfRecruitment.eps.....}{density panel plots of annual recuitment,} \cr
  \code{selectivityMCMC.eps....}{\bold{not currently implemented},} \cr
  \code{traceBiomass.eps.......}{panel plots of annual spawning biomass with running quantiles,} \cr
  \code{traceRecruits.eps......}{panel plots of annual recruitmwent traces with running quantiles,} \cr
  \code{traceParams.eps........}{panel plots of parameter traces with running quantiles,} \cr
  \code{splitChain.eps.........}{panel plots of cumulative parameter estimate chains,} \cr
  \code{VBcatch.eps............}{boxplots of annual vulnerable biomass and barplots of catch,} \cr
  \code{BVBnorm.eps............}{spawning and vulnerable biomass relative to their virgin levels,} \cr
  \code{Bproj.eps..............}{boxplots of spawning biomass -- MCMCs and projections,} \cr
  \code{Rproj.eps..............}{boxplots of recruitment -- MCMCs and projections,} \cr
  \code{RprojOnePolicy.eps.....}{boxplots of recruitment for one catch policy,} \cr
  \code{snail.eps..............}{time series of \eqn{u_t/u_{MSY}}{ut/umsy} \emph{vs.} \eqn{B_t/B_{MSY}}{Bt/Bmsy},} \cr
  \code{pairs[1,...,n].eps.....}{pairs plot of parameter MCMC samples (6 pars/page).} \cr
  \code{pairsMSY.eps...........}{pairs plot of MCMC samples for MSY-related parameters.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-25}
}
\note{
  The function is currently set up to send figures to encapsulated postscript (eps)
  and/or portable network graphics (png) files only. If a user wants other file types
  (e.g., windows metafiles), then additional lines of code need to be added to invoke
  the appropriate device(s).
}
\seealso{
  \code{\link{plt.mcmcGraphs}}, \code{\link{plotIndexNotLattice}}, 
  \code{\link{plotCPUE}}, \code{\link{plotChains}}, \code{\link{plotSnail}}
}
\keyword{hplot}
\keyword{connection}% __ONLY ONE__ keyword per line


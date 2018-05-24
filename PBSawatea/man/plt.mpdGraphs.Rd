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
   colGear=rep(c("black","blue"),ngear)[1:ngear])

plt.mcmcGraphs(mcmcObj, projObj=NULL, mpdObj=NULL, save=FALSE, 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, ngear=1,
   ylim.recruitsMCMC=NULL, ylim.exploitMCMC=NULL,
   ylim.VBcatch=NULL, ylim.BVBnorm=NULL,
   xlim.snail=NULL, ylim.snail=NULL,
   plotPolicies=names(projObj$Y[1:6]),
   onePolicy=names(projObj$Y[2]), mpd=list(),
   SAR.width=7.5, SAR.height=4, trevObj=NULL)
}
\arguments{
  \item{obj}{an Awatea MPD results object (\emph{e.g.}, \code{currentRes})}
  \item{save}{\bold{not used}}
  \item{ssnames}{survey series names for figure labelling and image creation}
  \item{ptypes}{plot types for the figures: \code{eps} = encapsulated postscript, \code{png} = portable network graphics.}
  \item{pngres}{resolution (pixels/inch) for the \code{png} file.}
  \item{ngear}{number of commercial gears}
  \item{pchGear, ltyGear, colGear}{symbols, line types, and colours for plotting data by commercial gear.}
  \item{mcmcObj}{an Awatea MCMC object (\emph{e.g.}, \code{currentMCMC})}
  \item{projObj}{an Awatea projected biomass object (\emph{e.g.}, \code{currentProj})}
  \item{mpdObj}{an Awatea MPD results object (\emph{e.g.}, \code{currentRes})}
  \item{ylim.recruitsMCMC}{y-limits for the recruitsMCMC plot}
  \item{ylim.exploitMCMC}{y-limits for the recruitsMCMC plot}
  \item{ylim.VBcatch}{y-limits for the VBcatch plot}
  \item{ylim.BVBnorm}{y-limits for the BVBnorm plot}
  \item{xlim.snail}{x-limits for the snail-trail plot}
  \item{ylim.snail}{y-limits for the snail-trail plot}
  \item{plotPolicies}{six-policies projections to plot}
  \item{onePolicy}{the one policy from the above six to use for some figures}
  \item{mpd}{list of MPD values for parameters, biomass, and recruitment}
  \item{SAR.width}{width (inches) of an output WMF file for SAR figures}
  \item{SAR.height}{height (inches) of an output WMF file for SAR figures}
  \item{trevObj}{various MCMC values suggested by Trevor Davies (Dalhousie) for comparison: \cr
    \code{c("h", "M_1", "M_2", "B0", "MSY", "Bmsy", "umsy")}. }
}
\details{
  Creates a whole heap of encapsulated postscript (eps) and portable network graphics (png) files.
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


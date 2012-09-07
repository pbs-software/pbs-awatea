\name{plt.mpdGraphs}
\alias{plt.mpdGraphs}
\alias{plt.mcmcGraphs}
\title{
  Plot a Set of Figures for MPD and MCMC
}
\description{
  Plot a set of postscript figures \code{.eps} for the MPD (mode of the posterior distribution)
  and MCMC (Monte Carlo Markoff Chain) results.
}
\usage{
plt.mpdGraphs(obj, save=FALSE, ssnames=paste("Ser",1:9,sep=""))
plt.mcmcGraphs(mcmcObj, projObj, save=FALSE, 
          ylim.recruitsMCMC = NULL, ylim.exploitMCMC = NULL,
          ylim.VBcatch = NULL, ylim.BVBnorm = NULL,
          xlim.snail = NULL, ylim.snail = NULL,
          plotPolicies = names(projObj$Y[1:6]),
          onePolicy=names(projObj$Y[2]), mpd = list())
}
\arguments{
  \item{obj}{an Awatea results object (\emph{e.g.}, \code{currentRes})}
  \item{mcmcObj}{an Awatea MCMC object (\emph{e.g.}, \code{currentMCMC})}
  \item{projObj}{an Awatea projected biomass object (\emph{e.g.}, \code{currentProj})}
  \item{save}{\bold{not used}}
  \item{ssnames}{survey series names for figure labelling and image creation}
  \item{ylim.recruitsMCMC}{y-limits for the recruitsMCMC plot}
  \item{ylim.exploitMCMC}{y-limits for the recruitsMCMC plot}
  \item{ylim.VBcatch}{y-limits for the VBcatch plot}
  \item{ylim.BVBnorm}{y-limits for the BVBnorm plot}
  \item{xlim.snail}{x-limits for the snail-trail plot}
  \item{ylim.snail}{y-limits for the snail-trail plot}
  \item{plotPolicies}{six-policies projections to plot}
  \item{onePolicy}{the one policy from the above six to use for some figures}
  \item{mpd}{list of MPD values for parameters, biomass, and recruitment.}
}
\details{
  Creates a whole heap o' postscript files.
}
\value{
  \bold{plt.mpdGraphs} creates the following postscript files: \cr
  \code{exploit.eps............}{annual exploitation rate,} \cr
  \code{recruits.eps...........}{annual recruitment at age 1,} \cr
  \code{selectivity.eps........}{selectivity curves for commercial gear(s) and survey(s),} \cr
  \code{ageComm.eps............}{fits to annual commercial age composition (panel plots),} \cr
  \code{ageSurv.eps............}{fits to annual survey age composition (panel plots),} \cr
  \code{survIndSer.eps.........}{four figures of survey indices (calls \code{\link{plotIndexNotLattice}}),} \cr
  \code{CPUEser.eps............}{CPUE indices with error bars (calls \code{\link{plotCPUE}}),} \cr
  \code{commAgeResids.eps......}{standardised residuals for commercial gear,} \cr
  \code{survAgeResidsSer.eps...}{standardised residuals for surveys,} \cr
  \code{meanAge.eps............}{mean age for catch and surveys,} \cr
  \code{stockRecruit.eps.......}{stock recruitment function.} \cr

  \bold{plt.mcmcGraphs} creates the following postscript files: \cr
  \code{recruitsMCMC.eps.......}{boxplots of annual MCMC recruitment,} \cr
  \code{exploitMCMC.eps........}{boxplots of annual MCMC exploitation rate,} \cr
  \code{pdfBiomass.eps.........}{density panel plots of annual female spawning biomass,} \cr
  \code{pdfRecruitment.eps.....}{density panel plots of annual recuitment,} \cr
  \code{selectivityMCMC.eps....}{\bold{not currently implemented},} \cr
  \code{traceRecruits.eps......}{panel plots of annual recruitmwent traces with running quantiles,} \cr
  \code{traceBiomass.eps.......}{panel plots of annual spawning biomass with running quantiles,} \cr
  \code{traceParams.eps........}{panel plots of parameter traces with running quantiles,} \cr
  \code{splitChain.eps.........}{panel plots of cumulative parameter estimate chains,} \cr
  \code{VBcatch.eps............}{boxplots of annual vulnerable biomass and barplots of catch,} \cr
  \code{BVBnorm.eps............}{spawning and vulnerable biomass relative to their virgin levels,} \cr
  \code{Bproj.eps..............}{boxplots of spawning biomass -- MCMCs and projections,} \cr
  \code{snail.eps..............}{time series of \eqn{u_t/u_{MSY}}{ut/umsy} \emph{vs.} \eqn{B_t/B_{MSY}}{Bt/Bmsy},} \cr
  \code{pairs[1,2,3].eps.......}{pairs plot of parameter MCMC samples.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  The function \code{plt.mpdGraphs} was recently cleaned up by AME. 
  Be wary not to use \code{currentRes} (globally available when the function 
  should be using \code{obj}, which is the results object passed to the function.

  The same is true for \code{plt.mcmcGraphs} where the global objects 
  \code{currentMCMC} and \code{currentProj} are use within the function
  rather than using the local objects \code{mcmcObj} and \code{projObj}.
}
\seealso{
  \code{\link{plt.mcmcGraphs}}, \code{\link{plotIndexNotLattice}}, 
  \code{\link{plotCPUE}}, \code{\link{plotChains}}, \code{\link{plotSnail}}
}
\keyword{hplot}
\keyword{connection}% __ONLY ONE__ keyword per line


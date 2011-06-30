\name{mainMenu}
\alias{mainMenu}
\alias{loadMenu}
\alias{mpdMenu}
\alias{mcmcMenu}
\alias{utilMenu}
\title{
  Create a Menu of Options and Actions
}
\description{
  From a main menu, the user can choose various options and actions. \cr
  We tend to do everything from the command line so the menu functionality
  fosters the warning: \emph{CAVEAT EMPTOR}.
}
\usage{
mainMenu()
loadMenu()
mpdMenu()
mcmcMenu()
utilMenu()
}
\details{
  \bold{Main menu items:} \cr
     \code{Import files} \cr 
     \code{MPD plots} \cr 
     \code{Plot all MPD graphs} \cr 
     \code{Save all MPD plots to PNG} \cr 
     \code{MCMC plots} \cr 
     \code{Plot all MCMC plots} \cr
     \code{Save all MCMC plots to PNG} \cr 
     \code{Close all graphics windows} \cr 
     \code{Help & Utilities}

  \bold{Load menu items:} \cr
     \code{Get Awatea res file} \cr 
     \code{Get Awatea MCMC file} \cr
     \code{Get Awatea projection file} \cr 
     \code{Load all res files in working directory} \cr
     \code{Get PJS Delay Difference MCMC+Projection} 

  \bold{MPD menu items:} \cr
    \code{Plot biomass, recruitment, catch} \cr 
    \code{Plot numbers at age} \cr
    \code{Plot selectivity and maturity} \cr 
    \code{Plot commercial catch-at-age results} \cr
    \code{Plot survey catch-at-age results} \cr 
    \code{Plot survey catch-at-length results} \cr
    \code{Plot abundance index} \cr 
    \code{All residual plots} \cr
    \code{Plot multi-panel biomass, recruitment, catch} \cr 
    \code{Plot multi-panel exploitation rate} \cr
    \code{Plot alternative numbers at age}

  \bold{MCMC menu items:} \cr
    \code{Plot biomass and projections by policy} \cr
    \code{Probability of projection biomass > reference} \cr
    \code{Expectation of projection biomass / reference} \cr
    \code{Plot biomass posterior densities (plotDens)} \cr
    \code{Plot recruitment posterior densities (plotDens)} \cr
    \code{Plot parameter posterior densities (plotDens)} \cr
    \code{Plot cumulative quantiles (plotCumu)} \cr
    \code{Plot traces (plotTrace)} \cr
    \code{Plot PJS traces (plt.allTraces)}

  \bold{Utils  menu items:} \cr
    \code{scape Help} \cr
    \code{scapeMCMC Help} \cr
    \code{Portrait graphsheet} \cr
    \code{Landscape graphsheet}
}
\seealso{
  \code{\link{get.resFile}}, \code{\link{importCol2}}
}
\keyword{interface}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

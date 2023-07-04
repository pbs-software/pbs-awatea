\name{plt.ageResidsPOP}
\alias{plt.ageResidsPOP}
\alias{plt.ageResidsqqPOP}
\alias{plt.cohortResids}
\alias{plt.yearResidsPOP}
\title{
  Plot Model Residuals
}
\description{
  Plot model residuals as quantile boxes.
}
\usage{
plt.ageResidsPOP(obj, ages=NULL, main=NULL, lang="e")

plt.ageResidsqqPOP(obj, ages=c(2,60),
   pct=c(5,25,50,75,95),  main=NULL, lang="e")

plt.yearResidsPOP (obj, ages=NULL, main=NULL,
   fill.in=TRUE, lang="e", ...)

plt.cohortResids(obj, ages=NULL, main=NULL, lang="e")
}
\arguments{
  \item{obj}{\code{data.frame} -- output from function \code{'stdRes.CA'}}
  \item{ages}{\code{numeric} -- age classes to plot}
  \item{pct}{\code{numeric} -- quantiles to show in quantile boxes or qq-plot}
  \item{main}{\code{character} -- title for plot if desired}
  \item{fill.in}{\code{logical} -- if \code{TRUE}, add missing years to quantile boxes}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french)}
  \item{...}{\code{dots} -- additional arguments for function \code{'quantbox'}.}
}
\details{
  \code{plt.ageResidsPOP.....}{plot age class residuals as quantile boxes.} \cr
  \code{plt.ageResidsqqPOP...}{plot age class residuals as qq-plot.} \cr
  \code{plt.yearResidsPOP....}{plot age residuals by year as quantile boxes.} \cr
  \code{plt.cohortResids.....}{plot age residuals by cohort as quantile boxes.} \cr
}
\author{
First created by:\cr
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew M. Edwards}, PBS, DFO, Nanaimo BC

Modified by:\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Regional Headquarters (RHQ), Vancouver BC\cr
  Last modified \code{Rd: 2023-06-28}
}
\note{
  Some trouble noted adding text and legend.
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{panelChains}}, 
  \code{\link[PBSawatea]{panelBoxes}}, 
  \code{\link[PBSawatea]{plotCPUE}},
  \code{\link[PBSawatea]{plotIndexNotLattice}},
  \code{\link[PBSawatea]{plt.allTraces}},
  \code{\link[PBSawatea]{plt.expRate}}, 
  \code{\link[PBSawatea]{plt.idx}},
  \code{\link[PBSawatea]{plt.mcmcGraphs}},
  \code{\link[PBSawatea]{stdRes.CA}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

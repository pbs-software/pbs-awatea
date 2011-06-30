\name{plt.ageResidsPOP}
\alias{plt.ageResidsPOP}
\alias{plt.ageResidsqqPOP}
\alias{plt.yearResidsPOP}
\alias{plt.cohortResids}
\title{
  Plot Model Residuals
}
\description{
  Plot model residuals as boxplots or qq-plots.
}
\usage{
plt.ageResidsPOP  (obj, ages=c(2,60), pct=c(5,25,50,75,95), main=NULL)
plt.ageResidsqqPOP(obj, ages=c(2,60), pct=c(5,25,50,75,95), main=NULL)
plt.yearResidsPOP (obj, ages=c(2,60), pct=c(5,25,50,75,95),
                   main=NULL, fill.in=TRUE, ... )
plt.cohortResids  (obj, ages=c(2,59), pct=c(5,25,50,75,95), main=NULL)

}
\arguments{
  \item{obj}{output from \code{\link{stdRes.CA}}.}
  \item{ages}{age classes to plot.}
  \item{pct}{quantiles to show in boxplot or qq-plot.}
  \item{main}{title for plot if desired.}
  \item{fill.in}{logical: if \code{TRUE}, add missing years to boxplot.}
  \item{...}{additional arguments for \code{boxplot}.}
}
\details{
  \code{plt.ageResidsPOP.....}{plot age class residuals as boxplots.} \cr
  \code{plt.ageResidsqqPOP...}{plot age class residuals as qq-plot.} \cr
  \code{plt.yearResidsPOP....}{plot age residuals by year as boxplots.} \cr
  \code{plt.cohortResids.....}{plot age residuals by cohort as boxplots.} \cr
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Some trouble noted adding text and legend.
}
\seealso{
  \code{\link{plt.allTraces}}, \code{\link{plt.expRate}}, 
  \code{\link{plt.idx}}, \code{\link{plt.mcmcGraphs}}, 
  \code{\link{plotIndexNotLattice}}, \code{\link{plotChains}}, 
  \code{\link{plotCPUE}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{importMCMC.ddiff}
\alias{importMCMC.ddiff}
\alias{importProj.ddiff}
\title{
  Import Functions for PJS Delay Difference Model
}
\description{
  Make a \pkg{plotMCMC} object identical in format to the result of 
  \code{importMCMC} (or \code{importProj}) from PJS delay 
  difference model output. \cr
  The difference is that \eqn{B} is biomass defined by a
  delay difference model.
}
\usage{
importMCMC.ddiff()
importProj.ddiff(yrVal="2006")
}
\arguments{
  \item{yrVal}{\code{character} -- year for delay-difference model}
}
\value{
  \code{importMCMC.ddiff} returns a list object containing: \cr
  \tabular{ll}{
    \code{L} \tab likelihood MCMCs\cr
    \code{P} \tab parameter MCMCs\cr
    \code{B} \tab spawning biomass MCMCs\cr
    \code{R} \tab recruitment MCMCs
  }

  \code{importProj.ddiff} returns a list object containing: \cr
  \tabular{ll}{
    \code{B} \tab projected biomass\cr
    \code{Y} \tab projected yields
  }
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC\cr
  Last modified \code{Rd: 2025-06-17}
}
\note{
  Get the biomass projection -- PJS does 1 year ahead projection.
  The column \code{"X"} appears as the last column because trailing 
  \code{","} exist in the \code{mcmcprojbiom.csv} file. \cr
  Note also that \code{"cat="} in \code{.csv} file becomes 
  \code{"cat."} in \code{read.table}.
}
\seealso{
In package \pkg{scape}:\cr
  \code{\link[scape]{importMCMC}},
  \code{\link[scape]{importProj}}\cr
In package \pkg{utils}:\cr
  \code{\link[utils]{read.table}}
}
% Programming
\keyword{file}
\keyword{interface}
\concept{PBSscape.r}


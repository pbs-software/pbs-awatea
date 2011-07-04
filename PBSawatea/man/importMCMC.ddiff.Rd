\name{importMCMC.ddiff}
\alias{importMCMC.ddiff}
\alias{importProj.ddiff}
\title{
  Import Functions for PJS Delay Difference Model
}
\description{
  Make a \pkg{scapeMCMC} object identical in format to the result of 
  \code{importMCMC} (or \code{importProj}) from PJS delay 
  difference model output. \cr
  The difference is that \eqn{B} is biomass defined by a
  delay difference model.
}
\usage{
importMCMC.ddiff()
importProj.ddiff()
}
\value{
  \code{importMCMC.ddiff} returns a list object containing: \cr
    \code{L...}{likelihood MCMCs,} \cr
    \code{P...}{parameter MCMCs,} \cr
    \code{B...}{spawning biomass MCMCs,} \cr
    \code{R...}{recruitment MCMCs.}

  \code{importProj.ddiff} returns a list object containing: \cr
    \code{B...}{projected biomass,} \cr
    \code{Y...}{projected yields.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\note{
  Get the biomass projection -- PJS does 1 year ahead projection.
  The column \code{"X"} appears as the last column because trailing 
  \code{","} exist in the \code{mcmcprojbiom.csv} file. \cr
  Note also that \code{"cat="} in \code{.csv} file becomes 
  \code{"cat."} in \code{read.table}.
}
\seealso{
  \code{\link[scapeMCMC]{importMCMC}}, \code{\link[scapeMCMC]{importProj}},
  \code{\link[utils]{read.table}}
}
% Programming
\keyword{file}
\keyword{interface}


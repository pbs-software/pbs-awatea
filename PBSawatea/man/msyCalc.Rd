\name{msyCalc}
\alias{msyCalc}
\title{
  Calculate the Maximum Sustainable Yield
}
\description{
  Load in \code{MSY.out} and calculate the MSY (maximum sustainable yield).
}
\usage{
msyCalc(dir=getwd(), error.rep=0, despike=FALSE)
}
\arguments{
  \item{dir}{\code{character} -- working directory}
  \item{error.rep}{\code{numeric} -- if \code{1}, report errors (reaching bounds), if \code{0} do not.}
  \item{despike}{\code{logical} -- if \code{TRUE}, attempt to remove false maxima (spikes) and use change in slope to locate MSY;
    if \code{FALSE}, use simple function \code{'max'} (original method).}
}
\value{
  Returns a list object containing: \cr
  \tabular{ll}{
    \code{yield} \tab maximum sustainable yield\cr
    \code{u}     \tab exploitation rate at MSY\cr
    \code{VB}    \tab vulnerable biomass at MSY\cr 
    \code{B}     \tab spawning biomass at MSY\cr
    \code{nProj} \tab number of projections needed to reach MSY
  }
}
\author{
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew Edwards}, Research Scientist, PBS, DFO, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-17}
}
\note{
  See \code{msyTestCreating.r} for full details when figuring this out.
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{findTarget}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{interation}
\keyword{manip}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{PBSscape.r}
% Use only one concept per line.
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

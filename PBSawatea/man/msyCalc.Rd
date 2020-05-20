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
    \code{yield...}{maximum sustainable yield,} \cr
    \code{u.......}{exploitation rate at MSY,} \cr
    \code{VB......}{vulnerable biomass at MSY,} \cr 
    \code{B.......}{spawning biomass at MSY,} \cr
    \code{nProj...}{numnber of projections needed to reach MSY.}
}
\author{
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew Edwards}, Research Scientist, PBS, DFO, Nanaimo BC\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2020-05-20}
}
\note{
  See \code{msyTestCreating.r} for full details when figuring this out.
}
\seealso{
  \code{\link{findTarget}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{interation}
\keyword{manip}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{MSY}
% Use only one concept per line.
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

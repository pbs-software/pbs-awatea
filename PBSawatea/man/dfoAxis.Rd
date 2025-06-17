\name{dfoAxis}
\alias{dfoAxis}
\title{
  Create DFO Axis
}
\description{
  Create a DFO axis that tries not to offend anyone, but offends everyone.
}
\usage{
dfoAxis(side, vec, prng=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{side}{\code{numeric} -- plot axis side (1, 2, 3, or 4)}
  \item{vec}{\code{numeric} -- vector of values that greatly exceed the axis range}
  \item{prng}{\code{numeric} -- appears to be a proportion outside the side's axis range used to clip the big vector}
}
\details{
  Looks suspiciously like the nonsense function \code{'dropLast'}
  but dressed up to impress someone.
  It actually makes more sense than \code{'dropLast'} alone in that ticks representing, say, 
  years 1900 to 2100 by 10 years will be clipped to an expansion of the axis' range of years.
}
\value{
  A DFO'd axis with less damage than a DFO'd computer.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\seealso{
\code{\link[PBSawatea]{dropLast}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{hogwash}
% \concept{ ~cpt2 }
% Use only one concept per line.

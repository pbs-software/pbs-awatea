\name{splitGear}
\alias{splitGear}
\title{
  Split Data by Gear
}
\description{
  Split data frames, where fields represent year-gear combinations, into one list object with gear-specific entries.
}
\usage{
splitGear(dat, fn=function(x) { sum(x, na.rm = TRUE) })
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{\code{numeric} -- usually a data frame of vulnerable biomass from MCMCs}
  \item{fn}{\code{function} -- not used for nada}
}
\details{
  For a data frame with fields that represent year-gear combinations, separate gear-specific fields and populate a list object by gear.
  Each list element comprises a data frame with dimensions of (MCMC sample number, year).
  Therefore, if there are three gear types in the original data frame, the list object will contain three gear-specific data frames.
  The function will also do this if there is only one gear type.
  Primarily used to separate gear-specific vulnerable biomass and exploitation rates from Awatea MCMC output.
}
\value{
  List object of gear-specific data frames.
}
\note{
  The function is not used much but it could be. 
  For instance, in \code{'plotSnail'} it could be used to split \code{'UoverUmsy'} (but it's not). C'est la vie.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Regional Headquarters (RHQ), Vancouver BC\cr
  Last modified \code{Rd: 2023-06-29}
}
\seealso{
\code{\link[PBSawatea]{plotTraj}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{MCMC}
% \concept{ ~cpt2 }
% Use only one concept per line.

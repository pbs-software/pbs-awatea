\name{repeatMPD}
\alias{repeatMPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Repeat MPDs for Likelihood Profiles
}
\description{
  Repeat MPDs to create likelihood profiles, usually for various levels (fixed) of \eqn{R_0}{R0} or \eqn{M}.
}
\usage{
repeatMPD(M=seq(0.07,0.08,0.01), A=c(40,45,50), R0=NULL, 
   prefix="WWR-MMM-", clean=FALSE, argsMPD="", dpY)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{\code{numeric} -- vector of natural mortality values}
  \item{A}{\code{numeric} -- vector of age plus-classes}
  \item{R0}{\code{numeric} -- vector of time-zero (equilibrium) recruitment values}
  \item{prefix}{\code{character} -- prefix for output names, usually species code
    and capitals designating the varying factor in the likelihood analysis}
  \item{clean}{\code{logical} -- if \code{TRUE}, delete the typical garbage produced by ADMB}
  \item{argsMPD}{\code{character} -- extra arguments to send to the ADMB call that fits the MPD}
  \item{dpY}{\code{numeric} -- number of decimal places to render the fixed factor as a unique character set.}
}
\details{
  Sets off a chain of MPD fits using an incremental factor.
}
\value{
  A series of results files (\code{.par}, \code{.res}, \code{.lik}, etc.)
  identified by species code and rendered factor level.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Regional Headquarters (RHQ), Vancouver BC\cr
  Last modified \code{Rd: 2023-06-29}
}
\seealso{
\code{\link[PBSawatea]{plotLikes}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{likelihood profiles}
% \concept{ ~cpt2 }
% Use only one concept per line.

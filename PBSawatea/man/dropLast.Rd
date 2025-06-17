\name{dropLast}
\alias{dropLast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Drop Last Value
}
\description{
  Drop the last value in a vector if it exceeds a target.
}
\usage{
dropLast(vec, target, prng=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{\code{numeric} -- vector of numbers (e.g., x-axis ticks)}
  \item{target}{\code{numeric} -- target number, maybe a limit (e.g., x-axis limit \code{par()$usr[2]})}
  \item{prng}{\code{numeric} -- proportion bracketing target}
}
\details{
  Not sure what this does, but it's used in PBSawatea for some reason.
}
\value{
  A clipped vector, perchance?
}
\note{
  Q'est-ce qu'il mange en hiver?
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Vivian Smith-Smythe-Smith}, Program Headless -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\seealso{
\code{\link[PBSawatea]{dfoAxis}}
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

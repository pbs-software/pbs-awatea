\name{plt.selectivity}
\alias{plt.selectivity}
\title{
  Plot Selectivity Fits
}
\description{
  Plot selectivity fits for each fleet (survey and/or commercial fishery)
}
\usage{
plt.selectivity(obj, mainTitle="Rockfish", ptypes=tcall(PBSawatea)$ptype,
   pngres=400, lang=c("e","f"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{\code{list} -- MPD object (usually \code{'currentRes'}) containing MLE fitted data}
  \item{mainTitle}{\code{character} -- title for selectivity plot (usually species name)}
  \item{ptypes}{\code{character} -- vector of plot types (e.g., \code{c("png","win")})}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) of the \code{png} file}
  \item{lang}{\code{character} -- vector of languages using first letter: \code{"e"}=english, \code{"f"}=french}
}
\details{
  Plots MLE fits for selectivity; female panels also show maturity ogive for reference.
}
\value{
  To die, to sleep -- to sleep, perchance to dream -- ay, there's the rub
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\seealso{
In package \pkg{PBSawatea}:\cr
\code{\link[PBSawatea]{plt.mpdGraphs}},
\code{\link[PBSawatea]{plotAges}}\cr
In package \pkg{PBSmodelling}:\cr
\code{\link[PBSmodelling]{tcall}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{hplot}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{plotFuns.r}
% \concept{ ~cpt2 }
% Use only one concept per line.

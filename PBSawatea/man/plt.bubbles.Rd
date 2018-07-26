\name{plt.bubbles}
\alias{plt.bubbles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Proportions-at-Age as Bubbles
}
\description{
  Plot proportions-at-age by representing annual age distributions as proportional bubbles
  for model inputs and fits by sex and gear or survey. 
}
\usage{
plt.bubbles(mpdObj, nsex=2, ptypes=tcall(PBSawatea)$ptype, 
   pngres=400, redo.Graphs=TRUE, lang=c("e","f"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpdObj}{\code{list} -- the MPD object created by the \code{importRes} function (e.g., \code{currentRes}).}
  \item{nsex}{\code{numeric} -- number of sexes in the model}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{"eps"} = encapsulated postscript, \code{"png"} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{.png} file.}
  \item{redo.Graphs}{\code{logical} -- if \code{TRUE}, re-create the figures.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
}
\details{
  The function \code{plt.bubbles} loops through all the catch@age (CA) 
  objects for commercial (CAc) and survey (CAs) observations and fits by 
  sex and gear/survey to produce bubble plots of the proportions-at age.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-26}
}
\seealso{
  \code{\link{plotCA}}, \code{\link{plt.recdev}}, \code{\link{plt.catch}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot} % __ONLY ONE__ keyword per line

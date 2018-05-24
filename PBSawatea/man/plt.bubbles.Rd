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
plt.bubbles(mpdObj, nsex=2, ptypes=tcall(PBSawatea)$ptype, pngres=400, redo.Graphs=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mpdObj}{the MPD object created by the \code{importRes} function (e.g., \code{currentRes}).}
  \item{nsex}{number of sexes in the model}
  \item{ptypes}{plot types for the figures: \code{eps} = encapsulated postscript, \code{png} = portable network graphics.}
  \item{pngres}{resolution (pixels/inch) for the \code{png} file.}
  \item{redo.Graphs}{\code{logical} -- if \code{TRUE}, re-create the figures.}
}
\details{
  The function \code{plt.bubbles} loops through all the catch@age (CA) 
  objects for commercial (CAc) and survey (CAs) observations and fits by 
  sex and gear/survey to produce bubble plots of the proportions-at age.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link{plotCA}}, \code{\link{plt.recdev}}, \code{\link{plt.catch}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot} % __ONLY ONE__ keyword per line

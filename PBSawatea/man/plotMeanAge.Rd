\name{plotMeanAge}
\alias{plotMeanAge}
\title{
  Plot Mean Age
}
\description{
  Plot model MPD estimates of mean age and observed mean age (with 95\% confidence interval).
}
\usage{
plotMeanAge(obj, useCA=TRUE, useSA=TRUE, CAnames, lang="e")
}
\arguments{
  \item{obj}{\code{list} -- The model's MPD object, usually \code{currentRes}.}
  \item{useCA}{\code{logical} -- if \code{TRUE}, then plot mean ages from the commercial fishery.}
  \item{useSA}{\code{logical} -- if \code{TRUE}, then plot mean ages from the surveys.}
  \item{CAnames}{\code{character} -- commercial fishery names, those that have composition (age) data.}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
}
\details{
  Provides a diagnostic plot, recommended by Francis (2011), that compares the model's
  calculated mean ages with observed mean ages for each abundance index.
}
\value{
  Contrary to CRAN's phobia about polluting the user's working environment, the function
  assigns a mean-age object \code{MA.pjs} to the global environment, \emph{and} dumps the object
  into an ASCII file (\code{MA.pjs.r}) and a binary file (\code{MA.pjs.rda}) in the MPD
  working directory. Stock assessments are messy and the user needs as much output as possible.
}
\references{
  Francis, R.I.C.C. (2011)
  Data weighting in statistical fisheries stock assessment models.
  \emph{Can. J. Fish. Aquat. Sci.} \bold{68}(6): 1124--1138.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\seealso{
\code{\link[PBSawatea]{importPar}},
\code{\link[PBSawatea]{importRes}},
\code{\link[PBSawatea]{MAfun}},
\code{\link[PBSawatea]{reweight}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")

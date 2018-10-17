\name{calc.refProbs}
\alias{calc.refProbs}
\alias{calc.refProbsHist}
\alias{calc.projExpect}
\alias{calc.projExpect2}
\alias{calc.projProbs}
\alias{calc.projProbs2}
\title{
  Calculate Reference Probabilities
}
\description{
  Calculate the expectation of projection to reference, and probability of 
  being greater than (or less than) the reference point.
}
\usage{
calc.refProbs(projObj=currentProj$B, refPlist=refPointsList, op=">")
calc.refProbsHist(projObj=currentProj$B, 
   refPlist=refPointsHistList[c("blimHRP","btarHRP")], op=">")
calc.projExpect (obj, projObj,  refYrs)
calc.projExpect2(obj, projObj, refList)
calc.projProbs  (obj, projObj,  refYrs)
calc.projProbs2 (obj, projObj, refList)
}
\arguments{
  \item{obj}{\code{matrix} -- biomass MCMCs.}
  \item{projObj}{\code{matrix} -- biomass projections.}
  \item{refYrs}{\code{numeric} -- vector of reference years}
  \item{refList}{\code{numeric list} -- reference years (numeric vectors).}
  \item{refPlist}{\code{numeric list} -- reference points.}
  \item{op}{\code{character} -- binary operator for comparing values in atomic vectors.}
}
\details{
  \code{calc.refProbs......}{Calculate the reference probabilities for decision tables.}\cr
  \code{calc.refProbsHist..}{Calculate the historical reference probabilities for decision tables.}\cr\cr
  \code{calc.projExpect....}{Calculate the expectation of projection to reference.}\cr
  \code{...................}{Compare reference years to projection years.}\cr
  \code{calc.projExpect2...}{Calculate expectation (projection biomass / reference biomass).}\cr\cr
  \code{calc.projProbs.....}{Calculate the probability of being greater than \code{refYrs}.}\cr
  \code{...................}{Compare reference years to projection years.}\cr
  \code{calc.projProbs2....}{Calculate the probability of being greater than \code{refYrs}.}\cr
  \code{...................}{Compare reference years to projection years.}
}
\value{
  Decision tables -- probability of projected biomass (or harvest rate)
  being greater than (or less than) a reference point.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2018-10-15}
}
\note{
  Functions probably deprecated because they are only called in \sQuote{\code{menu.r}}:\cr
  \code{calc.projExpect}, \code{calc.projExpect2}, \code{calc.projProbs}, \code{calc.projProbs2}
}
\seealso{
  \code{\link[PBSawatea]{calc.refVal}}
}
\keyword{array}
\keyword{list}


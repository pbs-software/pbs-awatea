\name{calc.projExpect}
\alias{calc.projExpect}
\alias{calc.projExpect2}
\alias{calc.projProbs}
\alias{calc.projProbs2}
\alias{calc.refProbs}
\alias{calc.refProbsHist}
\title{
  Calculate Expectations and Probabilities
}
\description{
  Calculate the expectation of projection to reference, and probability of being greater than reference.
}
\usage{
calc.projExpect  ( obj, projObj, refYrs  )
calc.projExpect2 ( obj, projObj, refList )
calc.projProbs   ( obj, projObj, refYrs  )
calc.projProbs2  ( obj, projObj, refList )
calc.refProbs    ( projObj=currentProj$B, refPlist=refPointsList )
calc.refProbsHist( projObj=currentProj$B, 
   refPlist=refPointsHistList[c("blimHRP","btarHRP")] )
}
\arguments{
  \item{obj}{matrix of biomass MCMCs.}
  \item{projObj}{matrix of biomass projections.}
  \item{refYrs}{numeric vector of reference years}
  \item{refList}{list of reference years (numeric vectors).}
  \item{refPlist}{list of reference points.}
}
\details{
  \code{calc.projExpect....}{Calculate the expectation of projection to reference.} \cr
  \code{...................}{Compare reference years to projection years.} \cr
  \code{calc.projExpect2...}{Calculate expectation (projection biomass / reference biomass).} \cr 
  \code{calc.projProbs.....}{Calculate the probability of being greater than \code{refYrs}.} \cr
  \code{...................}{Compare reference years to projection years.} \cr
  \code{calc.projProbs2....}{Calculate the probability of being greater than \code{refYrs}.} \cr
  \code{...................}{Compare reference years to projection years.} \cr
  \code{calc.refProbs......}{Calculate the reference probabilities (based on \code{calc.projProbs2}).}
  \code{calc.refProbsHist..}{Calculate the historical reference probabilities.}
}
\value{
  Decision tables
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link[PBSawatea]{calc.refVal}}
}
\keyword{array}
\keyword{list}


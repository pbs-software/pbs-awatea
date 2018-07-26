\name{plotSnail}
\alias{plotSnail}
\title{
  Phase Plot of Exploitation vs. Biomass
}
\description{
  Plot the historical progression of the ratio \eqn{u_t/u_{MSY}}{ut/umsy}
  against \eqn{B_t/B_{MSY}}{Bt/Bmsy}. Traditionally known as \sQuote{phase plots},
  some wine connoiseurs call these plots \emph{snail trails}.
}
\usage{
plotSnail(BoverBmsy, UoverUmsy, p=c(0.1,0.9), xLim=NULL, yLim=NULL,
   Lwd=1.5, ngear=1, assYrs=2011, outs=FALSE, lang="e")
}
\arguments{
  \item{BoverBmsy}{\code{numeric} -- matrix of \eqn{B_t}{Bt} over \eqn{B_{MSY}}{Bmsy}.}
  \item{UoverUmsy}{\code{numeric} -- matrix of \eqn{u_t}{ut} over \eqn{u_{MSY}}{umsy}.}
  \item{p}{\code{numeric} -- quantiles to show the bulk of the distribution.}
  \item{xLim}{\code{numeric} -- limits of the x-axis.}
  \item{yLim}{\code{numeric} -- limits of the y-axis.}
  \item{Lwd}{\code{numeric} -- line width of the snail trail.}
  \item{ngear}{\code{numeric} -- number of commercial gears}
  \item{assYrs}{\code{numeric} -- \sQuote{current} model years of past assessments.}
  \item{outs}{\code{logical} -- if \code{TRUE}, extend the range of the 
    quantile box whiskers to include the outliers (all data).}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french).}
}
\details{
  The graph attempts to show the time history of the exploitation rate 
  compared to the spawning biomass using a precautionary framework recast
  in Bayesian terms.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-26}
}
\note{
  The term \emph{snail trail} comes from PJS.
}
\seealso{
  \code{\link{plotBVBnorm}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

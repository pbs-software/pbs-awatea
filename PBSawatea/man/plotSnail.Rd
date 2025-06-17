\name{plotSnail}
\alias{plotSnail}
\title{
  Phase Plot of Exploitation vs. Biomass
}
\description{
  Plot the historical progression of the ratio \eqn{u_{t-1}/u_\mathrm{MSY}}{u[t-1]/u[MSY]}
  against \eqn{B_t/B_\mathrm{MSY}}{B[t]/B[MSY]}. These plots are known as 
  \sQuote{phase plots} by traditionalists, \sQuote{Kobe plots} by 
  tuna commissions, and \sQuote{snail trails} by NZ wine connoisseurs.
}
\usage{
plotSnail(BoverBmsy, UoverUmsy, p=c(0.05,0.95), xLim=NULL, yLim=NULL,
   Lwd=1.5, ngear=1, currYear=2021,  assYrs=NULL, outs=FALSE, Cnames, lang="e")
}
\arguments{
  \item{BoverBmsy}{\code{numeric} -- matrix of \eqn{B_t}{Bt} over \eqn{B_{MSY}}{Bmsy}}
  \item{UoverUmsy}{\code{numeric} -- matrix of \eqn{u_t}{ut} over \eqn{u_{MSY}}{umsy}}
  \item{p}{\code{numeric} -- quantiles to show the bulk of the distribution.}
  \item{xLim}{\code{numeric} -- limits of the x-axis}
  \item{yLim}{\code{numeric} -- limits of the y-axis}
  \item{Lwd}{\code{numeric} -- line width of the snail trail}
  \item{ngear}{\code{numeric} -- number of commercial gears}
  \item{currYear}{\code{numeric} -- current model year}
  \item{assYrs}{\code{numeric} -- model years of past assessments}
  \item{outs}{\code{logical} -- if \code{TRUE}, extend the range of the 
    quantile box whiskers to include the outliers (all data)}
  \item{Cnames}{\code{character} -- vector of commercial gear or fleet names}
  \item{lang}{\code{character} -- a letter that denotes the language for output: currently only \code{"e"} (english) or \code{"f"} (french)}
}
\details{
  The plot shows the time history of the exploitation rate compared
  to the spawning biomass using a precautionary framework recast in
  Bayesian terms.
}
\author{
First created by:\cr
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew M. Edwards}, PBS, DFO, Nanaimo BC

Modified by:\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\note{
  The term \emph{snail trail} comes from PJS.
  Younger cohorts have adopted the term \emph{Kobe plot}, so named by the 
  first joint meeting of the tuna regional fisheries management
  organizations (RFMOs), held in Kobe, Japan, in January 2007
  (Anonymous 2007)
}
\references{
  Anonymous (2007)
  \href{https://www.tuna-org.org/documents/other/kobe\%20report\%20english-appendices.pdf}{%
  Report of the Joint Meeting of Tuna RFMOs, January 20-26, 2007, Kobe, Japan},
  Participants: IATTC, ICCAT, IOTC, WCPFC, CCSBT.
}
\seealso{
\pkg{PBSawatea}:
  \code{\link[PBSawatea]{plotBVBnorm}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

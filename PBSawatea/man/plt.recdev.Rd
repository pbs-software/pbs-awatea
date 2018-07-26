\name{plt.recdev}
\alias{plt.recdev}
\alias{plt.recdevacf}
\alias{plt.initagedev}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Recruitment Deviations and ACF
}
\description{
  Plot recruitment deviations, the auto-correlation function of
  recuritment deviations, or the initial age deviations.
}
\usage{
plt.recdev(logRecDev, xint=5,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"))

plt.recdevacf(logRecDev, muC, logvC, A, years, yr1,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, redo.Graphs=TRUE, lang=c("e","f"))

plt.initagedev(logInitAgeDev, 
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logRecDev}{\code{numeric} -- vector of log recruitment deviations created 
    by the \code{importRes} function (e.g., \code{currentRes$Dev$Annual}).}
  \item{logInitAgeDev}{\code{numeric} -- vector of log initial age deviations created 
    by the \code{importRes} function (e.g., \code{currentRes$Dev$Initial}).}
  \item{muC}{\code{numeric} -- age at full selectivity for commercial gear created by
    the \code{importRes} function (e.g., \code{currentRes$extra$parameters$Sfullest}).}
  \item{logvC}{\code{numeric} -- log variance of commercial gear selectivity (left side) created by
    the \code{importRes} function (e.g., \code{currentRes$extra$parameters$log_varLest}).}
  \item{A}{\code{numeric} -- maximum age in the model (plus class)}
  \item{years}{\code{numeric} -- years covered by the population model reconstruction}
  \item{yr1}{\code{numeric} -- first year of age data}
  \item{xint}{\code{numeric} -- x-axis (year) intervals for small, unlabelled ticks}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{"eps"} = encapsulated postscript, \code{"png"} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{.png} file.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
  \item{redo.Graphs}{\code{logical} -- if \code{TRUE}, re-create the figures.}
}
\details{
  The function \code{plt.recdev} plots annual log recuritment deviations.

  The function \code{plt.recdevacf} plots the auto-correlation function (first 30 lags) for the log recuritment deviations.

  The function \code{plt.initagedev} plots the initial age deviation for each age in the model.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-26}
}
\seealso{
  \code{\link{plt.catch}}, \code{\link{plt.cpue}}, \code{\link{plt.bubbles}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot} % __ONLY ONE__ keyword per line

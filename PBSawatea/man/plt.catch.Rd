\name{plt.catch}
\alias{plt.catch}
\alias{plt.biomass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Catch and Biomass History
}
\description{
  Plot catch and/or model biomass by year. 
}
\usage{
plt.catch(years, Ct, xint=5, yint=250,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"))

plt.biomass(years, Bt, xint=5, yint=2500,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, 
   pname="Bt", xlab="Year", ylab="Spawning biomass (t), Bt", lang=c("e","f"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{years}{\code{numeric} -- vector of years}
  \item{Ct}{\code{numeric} -- vector or matrix of annual catches; if a matrix, then rows=years and cols=gear.}
  \item{Bt}{\code{numeric} -- vector or matrix of annual biomass; if a matrix, then rows=years and cols=gear.}
  \item{xint}{\code{numeric} -- x-axis (year) intervals for small, unlabelled ticks}
  \item{yint}{\code{numeric} -- y-axis (catch or biomass) intervals for small, unlabelled ticks}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{eps} = encapsulated postscript, \code{png} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{png} file.}
  \item{pname}{\code{character} -- prefix name for an image file}
  \item{xlab}{\code{character} -- label for the x-axis}
  \item{ylab}{\code{character} -- label for the y-axis}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
}
\details{
  The function \code{plt.catch} plots catch as vertical bars and spits out 
  two files named \code{catch.eps} and \code{catchsmall.eps} (and/or \code{.png} files).

  The function \code{plt.biomass} plots biomass as points and allows more flexibility to handle 
  various biomass inputs (e.g., \eqn{B_t}, \eqn{B_t/B_0}, \eqn{V_{tg}/V_{0g}}).
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\seealso{
  \code{\link{plt.cpue}}, \code{\link{plt.recdev}}, \code{\link{plt.bubbles}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot} % __ONLY ONE__ keyword per line

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
   ptypes=tcall(PBSawatea)$ptype, pngres=400)

plt.biomass(years, Bt, xint=5, yint=2500,
   ptypes=tcall(PBSawatea)$ptype, pngres=400, 
   pname="Bt", xlab="Year", ylab="Spawning biomass (t), Bt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{years}{numeric vector of years}
  \item{Ct}{vector or matrix of annual catches; if a matrix, then rows=years and cols=gear.}
  \item{Bt}{vector or matrix of annual biomass; if a matrix, then rows=years and cols=gear.}
  \item{xint}{x-axis (year) intervals for small, unlabelled ticks}
  \item{yint}{y-axis (catch or biomass) intervals for small, unlabelled ticks}
  \item{ptypes}{plot types for the figures: \code{eps} = encapsulated postscript, \code{png} = portable network graphics.}
  \item{pngres}{resolution (pixels/inch) for the \code{png} file.}
  \item{pname}{prefix name for an image file}
  \item{xlab}{label for the x-axis}
  \item{ylab}{label for the y-axis}
}
\details{
  The function \code{plt.catch} plots catch as vertical bars and spits out 
  two files named \code{catch.eps} and \code{catchsmall.eps} (and/or \code{.png} files).

  The function \code{plt.biomass} plots biomass as points and allows more flexibility to handle 
  various biomass inputs (e.g., \eqn{B_t}, \eqn{B_t/B_0}, \eqn{V_{tg}/V_{0g}}).
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link{plt.cpue}}, \code{\link{plt.recdev}}, \code{\link{plt.bubbles}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot} % __ONLY ONE__ keyword per line

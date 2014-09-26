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
   ptypes = c("eps", "png"), pngres = 150)

plt.recdevacf(logRecDev, muC, logvC, A, years, yr1,
   ptypes=c("eps","png"), pngres=150)

plt.initagedev(logInitAgeDev, 
   ptypes=c("eps","png"), pngres=150 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{logRecDev}{vector of log recruitment deviations created 
    by the \code{importRes} function (e.g., \code{currentRes$Dev$Annual}).}
  \item{logInitAgeDev}{vector of log initial age deviations created 
    by the \code{importRes} function (e.g., \code{currentRes$Dev$Initial}).}
  \item{muC}{age at full selectivity for commercial gear created by
    the \code{importRes} function (e.g., \code{currentRes$extra$parameters$Sfullest}).}
  \item{logvC}{log variance of commercial gear selectivity (left side) created by
    the \code{importRes} function (e.g., \code{currentRes$extra$parameters$log_varLest}).}
  \item{A}{maximum age in the model (plus class)}
  \item{years}{years covered by the population model reconstruction}
  \item{yr1}{first year of age data}
  \item{xint}{x-axis (year) intervals for small, unlabelled ticks}
  \item{ptypes}{plot types for the figures: \code{eps} = encapsulated postscript, \code{png} = portable network graphics.}
  \item{pngres}{resolution (pixels/inch) for the \code{png} file.}
}
\details{
  The function \code{plt.recdev} plots annual log recuritment deviations.

  The function \code{plt.recdevacf} plots the auto-correlation function (first 30 lags) for the log recuritment deviations.

  The function \code{plt.initagedev} plots the initial age deviation for each age in the model.
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link{plt.catch}}, \code{\link{plt.cpue}}, \code{\link{plt.bubbles}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot} % __ONLY ONE__ keyword per line

\name{mochaLatte}
\alias{mochaLatte}
\title{
  Mock Lattice with a dash of nutmeg
}
\description{
  An attempt to simulate a \code{lattice} panel plot and gain greater control over
  the look of the figure without having to tweak endless sub-options.
}
\usage{
mochaLatte(dat, xfld, yfld, ffld, panel,
   strip=list(col=lucent("black",0.5), bg=lucent("moccasin",0.5), height=0.1, cex=1.4), 
   fn.ylim=range, ...)
}
\arguments{
  \item{dat}{\code{data.frame} -- object with fields that can be used as factors (without being factors) and fields with data to plot.}
  \item{xfld}{\code{character} -- name of the field with x-values.}
  \item{yfld}{\code{character} -- name of the field with y-values.}
  \item{ffld}{\code{character} -- name of the field with \sQuote{factors} that will be represented in each panel.}
  \item{panel}{\code{function} -- panel function that plots the data in some fashion.}
  \item{strip}{\code{list} -- controls for strip along the top of each panel: 
    outline colour (\code{col}), background colour (\code{bg}), height, and font size (\code{cex}).}
  \item{fn.ylim}{\code{function} -- to determine the y-limits of the panel plots (default function = \code{'range'}).}
  \item{\dots}{\code{list} -- aditional arguments for \code{plot}, \code{axis}, and \code{mtext}.}
}
\details{
  The function does not mock \code{lattice}, rather, it provides a \emph{faux} 
  \code{lattice}-like, panel-plotting framework that will undoubtedly evolve as it is used
  to replace other calls to \code{lattice} plots. There is no guarantee that this
  function is universally functional. But perhaps more importantly, it is an 
  alternative to embracing \code{ggplot}.
  Ain't nobody got time for that.
}
\value{
  Currently returns nothing but joy.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2019-05-31}
}
\note{
  The \code{lattice} package provides a versatile system that has evolved over many years;
  however, there is only so much time in a day that one can fiddle about with lattice options.
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{panelChains}},
  \code{\link[PBSawatea]{plotDensPOP}},
  \code{\link[PBSawatea]{plotIndexNotLattice}},
  \code{\link[PBSawatea]{plotTracePOP}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")

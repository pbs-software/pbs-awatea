\name{plt.cpue}
\alias{plt.cpue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot CPUE Fits
}
\description{
  Plot CPUE fits to annual indices (only if the model uses CPUE). 
}
\usage{
plt.cpue(cpueObj, ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cpueObj}{\code{data.frame} -- a CPUE object created by the \code{importRes} function (e.g., \code{currentRes$CPUE}).}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{"eps"} = encapsulated postscript, \code{"png"} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{.png} file.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
}
\details{
  The function \code{plt.cpue} is only activated in the Sweave routine if the
  Awatea model uses CPUE.
  The image file created is called \code{CPUEfit.eps} (and/or \code{.png}) and 
  shows the model fit through the CPUE points.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\seealso{
  \code{\link{plt.catch}}, \code{\link{plt.recdev}}, \code{\link{plt.bubbles}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot} % __ONLY ONE__ keyword per line

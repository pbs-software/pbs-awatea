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
plt.cpue(cpueObj, ptypes=tcall(PBSawatea)$ptype, pngres=400)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cpueObj}{a CPUE object created by the \code{importRes} function (e.g., \code{currentRes$CPUE}).}
  \item{ptypes}{plot types for the figures: \code{eps} = encapsulated postscript, \code{png} = portable network graphics.}
  \item{pngres}{resolution (pixels/inch) for the \code{png} file.}
}
\details{
  The function \code{plt.cpue} is only activated in the Sweave routine if the
  Awatea model uses CPUE.
  The image file created is called \code{CPUEfit.eps} (and/or \code{.png}) and 
  shows the model fit through the CPUE points.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link{plt.catch}}, \code{\link{plt.recdev}}, \code{\link{plt.bubbles}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot} % __ONLY ONE__ keyword per line

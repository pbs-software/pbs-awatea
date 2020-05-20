\name{redo.currentProj}
\alias{redo.currentProj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Recalculate the currentProj Object(s)
}
\description{
  Recalculate objects \code{'currentProj'} (CC policy) and \code{'currentProj2'} (HR policy).
}
\usage{
redo.currentProj(prj.dir, Ngear, sigmaR=0.9, mcsub=201:1200,
   compile.only=FALSE, recalc.proj=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prj.dir}{\code{character} -- the projection directory, usually under an MCMC directory and prefixed by \code{'PRJ'}.}
  \item{Ngear}{\code{numeric} -- number of gear types used in the model.}
  \item{sigmaR}{\code{numeric} -- standard deviation of recruitment residuals used in the model.}
  \item{mcsub}{\code{numeric} -- vector of samples used for MCMC
    (e.g., \code{201:1200} for \code{1000} samples, discarding the first \code{200} as burn-in).}
  \item{compile.only}{\code{logical} -- if \code{TRUE}, only evaluate the ADMB expression in the relevant projection directory:
    (e.g., \code{'awatea -ind BSR-2F.28.01.txt -mceval'}. }
  \item{recalc.proj}{\code{logical} -- if \code{TRUE}, recalculate the object(s)
    \code{'currentProj'} (for a CC policy) and/or \code{'currentProj2'} (for an HR policy).}
}
\details{
  Essentially repeats code that was use in \code{'run-SweaveMCMC.Snw'} so that user can
  recalculate projected catch policies if they were changed from the original.
}
\value{
	Objects \code{'currentProj'} (for a CC policy) and/or \code{'currentProj2'} (for an HR policy) in the MCMC directory.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2020-05-20}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{runSweaveMCMC}},
  \code{\link[PBSawatea]{importProjRec}},
  \code{\link[PBSawatea]{srFun}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{manip}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{Projection}
% \concept{ ~cpt2 }
% Use only one concept per line.

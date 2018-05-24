\name{plotACFs}
\alias{plotACFs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot Autocorrelations for MCMC Parameters
}
\description{
  Plot ACFs for a set of estimated parameters from MCMC samples.
}
\usage{
plotACFs(mcmcObj, lag.max = 60)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mcmcObj}{\code{list} -- collection of data frames from the MCMC, including \code{P} used by this function.}
  \item{lag.max}{\code{numeric} -- maximum lag at which to calculate the autocorelation function (see \code{\link[stats]{acf}}).}
}
\details{
  Extracts the estimated parameters from the data frame \code{currentMCMC$P} and 
  displays the first \code{lag.max} lags using the \code{acf} function from the \code{stats} package.
}
\value{
  Returns nothing at present.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-05-24}
}
\note{
  \emph{Caution} -- (from \code{plot.acf}) The confidence interval plotted in \code{plot.acf}
  is based on an \emph{uncorrelated} series and should be treated with appropriate caution. 
  Using \code{ci.type = "ma"} may be less potentially misleading. 
}
\seealso{
\code{\link[PBSawatea]{plotChains}},
\code{\link[PBSawatea]{plotCI}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")

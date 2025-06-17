\name{MAfun}
\alias{MAfun}
\title{
  Calculate Mean Age by Year
}
\description{
  Calculate mean ages from proportions-at-age 
  (modified from a subfunction in \code{runADMB}).
}
\usage{
MAfun(padata, brks=NULL)
}
\arguments{
  \item{padata}{\code{numeric} -- proportion-at-age data \code{CAc} or
    \code{CAs} from a call to \code{importRes}}
  \item{brks}{\code{numeric} -- breaks specified as numeric years to
    split the commercial data up into regimes that may account for
    index discontinuities (\bold{not used})}
}
\details{
  Mean age function supplied by Chris Francis (2011). \cr
  \code{padata} has fields: \cr
  \tabular{ll}{
    \code{Series} \tab series identifier\cr
    \code{Year}   \tab numeric year\cr
    \code{Age}    \tab age bin\cr
    \code{Obs}    \tab observed proportions\cr
    \code{Fit}    \tab predicted (fitted) proportions\cr
    \code{SS}     \tab sample size (effective \eqn{N})
  }
}
\value{
  List object of observed and expected mean ages, variance of expected ages,
  and a few bits and bobs.
}
\references{
Francis, R.I.C.C. (2011) 
  Data weighting in statistical fisheries stock assessment models. 
  \emph{Canadian Journal of Fisheries and Aquatic Sciences} \bold{68}: 1124-1138.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-17}
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{runADMB}},
  \code{\link[PBSawatea]{importRes}}
}
\keyword{univar}
\keyword{list}
\concept{utilFuns.r}


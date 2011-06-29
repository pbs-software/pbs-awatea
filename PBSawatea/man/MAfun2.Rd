\name{MAfun2}
\alias{MAfun2}
\title{
  Calculate Mean Age by Year
}
\description{
  Calculate mean ages from proportions-at-age 
  (modified from a subfunction in \code{runADMB}).
}
\usage{
MAfun2(padata, brks=NULL)
}
\arguments{
  \item{padata}{proportion-at-age data \code{CAc} or \code{CAs} from a 
    call to \code{importCol2}.}
  \item{brks}{breaks specified as numeric years to split the 
    commercial data up into regimes that may account for index 
    discontinuities (\bold{not used}). }
}
\details{
  Mean age function supplied by Chris Francis (2011). \cr
  \code{padata} has fields: \cr
  \code{Series...}{series identifier} \cr
  \code{Year.....}{numeric year} \cr
  \code{Age......}{age bin} \cr
  \code{Obs......}{observed proportions} \cr
  \code{Fit......}{predicted (fitted) proportions} \cr
  \code{SS.......}{sample size (effective \eqn{N})}
}
\value{
  List object of observed and expected mean ages, variance of expected ages,
  and a few bits and bobs.
}
\references{
  Francis, R.I.C.C. (2011, in press)
  Data weighting in statistical fisheries stock assessment models. 
  \emph{Canadian Journal of Fisheries and Aquatic Sciences}.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC.
}
\seealso{
  \code{\link{runADMB}}, \code{\link{importCol2}}
}
\keyword{univar}
\keyword{list}


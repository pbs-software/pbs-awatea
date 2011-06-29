\name{reweight}
\alias{reweight}
\title{
  Reweight Abundance and Composition Data
}
\description{
  Weight the abundance data by adjusting survey and CPUE CVs, and 
  weight the composition data by adjusting the effective sample size \eqn{N}.
}
\usage{
reweight(obj, cvpro=FALSE, mean.age=TRUE, ...)
}
\arguments{
  \item{obj}{an \code{AWATEA} class object created initially by \code{readAD}.}
  \item{cvpro}{CV process error added to CV observation error: \eqn{c_t = \sqrt{c_o^2 + c_p^2}}{c[t] = sqrt(c[o]^2 + c[p]^2)}; \cr
    if \code{FALSE} index CVs are reweighted using the standard deviation of normalized residuals.}
  \item{mean.age}{logical: if \code{TRUE}, use mean-age residuals to reweight
    the effective \eqn{N} for the age composition data (see Francis 2011); \cr
    if \code{FALSE}, reweight \eqn{N} using \eqn{ \Sigma (P(1-P)) / \Sigma (O-P)^2 }, 
    where \eqn{O} = observed proportions-at age and \eqn{P} = predicted/fitted proportions-at-age.}
  \item{\dots}{additional arguments to \code{reweight}.}
}
\details{
  For the reweight to work, a corresponding Awatea results file (\code{.res}) with the 
  same prefix as the input file must be available in the working directory before calling 
  \code{readAD}. This will populate the \code{output} slot with fitted data that the reweight needs.
}
\value{
  An \code{AWATEA} class cobject with the slots outlined in \code{\link{readAD}} 
  with the following slot populated by this function: \cr
  \code{reweight...}{list of reweight results:} \cr
  \code{..nrwt.....}{the number of the current reweighting,} \cr
  \code{..survey...}{survey indices with CV values (observed, fitted, normalised residuals, reweighted),} \cr
  \code{..cpue.....}{if used in the model, CPUE indices with CV values (as above),} \cr
  \code{..wNcpa....}{reweigthed effective \eqn{N} for commercial compositions (proportions-at-age),} \cr
  \code{..wNspa....}{reweigthed effective \eqn{N} for survey compositions,} \cr
  \code{..SDNR.....}{standard deviation of normalised residuals for abundance and composition data,} \cr
  \code{..wj.......}{weights for composition data from a mean-age weighted calculation (Francis 2011).}
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
  \code{\link{runADMB}}, \code{\link{readAD}}, \code{\link{stdRes.CA}}, \code{\link{stdRes.index}}
}
\keyword{models}
\keyword{methods}


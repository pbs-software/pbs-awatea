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
reweight(obj, A.rwt=1, C.rwt=1, cvpro=0, ...)
}
\arguments{
  \item{obj}{\code{S4 object} -- an \code{AWATEA} class object created initially by \code{'readAD'}.}
  \item{A.rwt}{\code{numeric} -- abundance reweight method: 0=no reweight, 1=add pocess error, 2=SDNR.}
  \item{C.rwt}{\code{numeric} -- composition reweight method: 0=no reweight, 1=Francis (2011) mean age, 2=SDNR.}
  \item{cvpro}{\code{numeric} -- CV process error added to CV observation error (see \bold{Details}); \cr
    if \code{cvpro} vector is shorter than \code{Nsurvey + Ncpue}, 
    it is repeated to match the number of survey series and CPUE series used to index abundance.}
  \item{\dots}{additional arguments to \code{reweight}.}
}
\details{
  For the reweight to work, a corresponding Awatea results file (\code{.res}) with the 
  same prefix as the input file must be available in the working directory before calling 
  \code{readAD}. This will populate the \code{output} slot with fitted data that the reweight needs.

  Reweighting of abundance data and composition data can now be performed independently of
  each other using the vectors \code{A.rwt} and \code{C.rwt}, respectively.
  This flexibility also means that a user could choose different reweight techniques
  for every reweight, but mixing techniques is not generally done.
  (See \code{\link[PBSawatea]{runADMB}} for more detail.)

  CV process error \eqn{c_\mathrm{p}}{c[p]} is typically added to observed error \eqn{c_\mathrm{o}}{c[o]} using:
  \deqn{c_\mathrm{t} = \sqrt{c_\mathrm{o}^2 + c_\mathrm{p}^2}}{c[t] = sqrt(c[o]^2 + c[p]^2)}
  For a method of calculating commercial CPUE \eqn{c_\mathrm{p}}{c[p]}, see \code{\link[PBSawatea]{splineCPUE}}.

  Composition reweighting adjusts the effective sample size (SS) of age frequency data.
  Generally, the mean age method of Francis (2011) is used.
  If \code{C.rwt=2}, SS are reweighted using the standard deviation of normalised residuals (SDNR):
  \deqn{\frac{\sum_{a=1}^A P_a(1-P_a)}{\sum_{a=1}^A (O_a-P_a)^2 }}{\Sigma[a=1]^A (P[a](1-P[a])) / \Sigma[a=1]^A (O[a]-P[a])^2 }
  where \eqn{a} = age, \eqn{A} = accumulator age class, \eqn{O_a}{O[a]} = observed proportions-at age and \eqn{P_a}{P[a]} = predicted/fitted proportions-at-age.
}
\value{
  An \code{AWATEA} class object with the slots outlined in \code{\link{readAD}} 
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
  Francis, R.I.C.C. (2011)
  \href{https://doi.org/10.1139/f2011-025}{Data weighting in statistical fisheries stock assessment models}.
  Can. J. Fish. Aquat. Sci. 68(6): 1124-1138.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Regional Headquarters (RHQ), Vancouver BC\cr
  Last modified \code{Rd: 2023-06-28}
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{readAD}},
  \code{\link[PBSawatea]{runADMB}},
  \code{\link[PBSawatea]{splineCPUE}},
  \code{\link[PBSawatea]{stdRes.CA}},
  \code{\link[PBSawatea]{stdRes.index}}
}
\keyword{models}
\keyword{methods}


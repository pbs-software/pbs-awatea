\name{runADMB}
\alias{runADMB}
\title{
  Run AD Model Builder Code for Awatea
}
\description{
  Run compiled AD Model Builder code called \code{Awatea} to reconstruct a
  population trajectory for a marine fish stock.
}
\usage{
runADMB(filename.ext, wd=getwd(), strSpp="XYZ", runNo=1, rwtNo=0,
     doMPD=FALSE, N.reweight=0, cvpro=FALSE, mean.age=TRUE, 
     doMCMC=FALSE, mcmc=1e6, mcsave=1e3, ADargs=NULL, verbose=FALSE, 
     doMSY=FALSE, msyMaxIter=15000., msyTolConv=0.01, endStrat=0.301, stepStrat=0.001,
     delim="-", awateaPath="E:/Projects/ADMB/Coleraine", ...)
}
\arguments{
  \item{filename.ext}{character file name including its extension.}
  \item{wd}{character string specifying the working directory, if different from the current working directory.}
  \item{strSpp}{string 3-letter code of the species.}
  \item{runNo}{the model run number.}
  \item{rwtNo}{the reweight number.}
  \item{doMPD}{logical: if \code{TRUE}, perform an MPD analysis.}
  \item{N.reweight}{the number of reweights to perform in the MPD analysis.}
  \item{cvpro}{CV process error added to CV observation error: \eqn{c_t = \sqrt{c_o^2 + c_p^2}}; \cr
    if \code{FALSE} index CVs are reweighted using the standard deviation of normalized residuals.}
  \item{mean.age}{logical: if \code{TRUE}, use mean-age residuals to reweight
    the effective \eqn{N} for the age composition data (see Francis 2011); \cr
    if \code{FALSE}, reweight \eqn{N} using \eqn{ \Sigma (P(1-P)) / \Sigma (O-P)^2 }, 
    where \eqn{O} = observed proportions-at age and \eqn{P} = predicted/fitted proportions-at-age.}
  \item{doMCMC}{logical: if \code{TRUE}, perform an MCMC analysis.}
  \item{mcmc}{number of MCMC iterations to perform.}
  \item{mcsave}{frequency of MCMC iterations to save.}
  \item{ADargs}{additional arguments for a call to Awatea.}
  \item{verbose}{logical: if \code{TRUE}, spew Awatea messages to the R console.}
  \item{doMSY}{logical: if \code{TRUE}, perform an MSY analysis.}
  \item{msyMaxIter}{maximum iterations for the MSY calculations.}
  \item{msyTolConv}{tolerance for convergence in the MSY calculations.}
  \item{endStrat}{maximum fishing mortality for the MSY analysis.}
  \item{stepStrat}{fishing mortality step size for the MSY analysis.}
  \item{delim}{character that delimits the components of the \code{filename} prefix.}
  \item{awateaPath}{system path on which the executable file \code{Awatea.exe} exists.}
  \item{\dots}{additional arguments (not currently used for any purpose).}
}
\details{
  This function is primarily used to automate MPD reweightings and to perform 
  the MSY calculations. The MCMCs are better run from a command line console on 
  a super computer as they generally require > 12 h to complete.
  
  Note that once an MCMC has been created, the user can run various projections separately using \cr 
  \code{awatea -ind filename.ext -mceval} \cr on the command line.
}
\references{
  Edwards, A.M., Starr, P.J., and Haigh, R. (2010, in revision). 
  Stock assessment for Pacific Ocean Perch (\emph{Sebastes alutus}) in Queen Charlotte Sound, British Columbia.
  \emph{Canadian Science Advisory Secretariat, Research Document}.

  Edwards, A.M., Haigh, R., and Starr, P.J. (2011, in revision). 
  Stock assessment and recovery potential assessment for Yellowmouth Rockfish (\emph{Sebastes reedi}) along the Pacific coast of Canada.
  \emph{Canadian Science Advisory Secretariat, Research Document}.
  
  Francis, R.I.C.C. (2011, in press)
  Data weighting in statistical fisheries stock assessment models. 
  \emph{Canadian Journal of Fisheries and Aquatic Sciences}.

  Hilborn, R., Maunder, M., Parma, A., Ernst, B. Payne, J., and Starr, P. (2003)
  Coleraine: a generalized age-structured stock assessment model. 
  School of Aquatic and Fishery Sciences, University of Washington, 54 p.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC.
}
\seealso{
  \code{\link{importRes}}, \code{\link{readAD}}, \code{\link{reweight}}, \code{\link{runSweave}}
}
\keyword{models}
\keyword{iteration}


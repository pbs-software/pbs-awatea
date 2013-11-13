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
runADMB(
  filename.ext, wd=getwd(), 
  strSpp="XYZ", runNo=1, rwtNo=0,
  doMPD=FALSE, N.reweight=0, cvpro=FALSE, mean.age=TRUE, 
  doMCMC=FALSE, mcmc=1e6, mcsave=1e3, ADargs=NULL, verbose=FALSE, 
  doMSY=FALSE, msyMaxIter=15000., msyTolConv=0.01, endStrat=0.301, stepStrat=0.001,
  delim="-", clean=FALSE, locode=FALSE, 
  awateaPath="C:/Users/haighr/Files/Projects/ADMB/Coleraine",
  codePath="C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
  ...
)
}
\arguments{
  \item{filename.ext}{\emph{character} file name including its extension.}
  \item{wd}{\emph{character} string specifying the working directory, if different from the current working directory.}
  \item{strSpp}{\emph{character} string 3-letter code of the species.}
  \item{runNo}{\emph{numeric}: the model run number.}
  \item{rwtNo}{\emph{numeric}: the reweight number.}
  \item{doMPD}{\emph{logical}: if \code{TRUE}, perform an MPD analysis.}
  \item{N.reweight}{\emph{numeric}: the number of reweights to perform in the MPD analysis.}
  \item{cvpro}{\emph{numeric/logical}: CV process error added to CV observation error: \eqn{c_t = \sqrt{c_o^2 + c_p^2}}; \cr
    if \code{FALSE} index CVs are reweighted using the standard deviation of normalized residuals.}
  \item{mean.age}{\emph{logical}: if \code{TRUE}, use mean-age residuals to reweight
    the effective \eqn{N} for the age composition data (see Francis 2011); \cr
    if \code{FALSE}, reweight \eqn{N} using \eqn{ \Sigma (P(1-P)) / \Sigma (O-P)^2 }, 
    where \eqn{O} = observed proportions-at age and \eqn{P} = predicted/fitted proportions-at-age.}
  \item{doMCMC}{\emph{logical}: if \code{TRUE}, perform an MCMC analysis.}
  \item{mcmc}{\emph{numeric}: number of MCMC iterations to perform.}
  \item{mcsave}{\emph{numeric}: frequency of MCMC iterations to save.}
  \item{ADargs}{\emph{character}: additional arguments for a call to Awatea.}
  \item{verbose}{\emph{logical}: if \code{TRUE}, spew Awatea messages to the R console.}
  \item{doMSY}{\emph{logical}: if \code{TRUE}, perform an MSY analysis.}
  \item{msyMaxIter}{\emph{numeric}: maximum iterations for the MSY calculations.}
  \item{msyTolConv}{\emph{numeric}: tolerance for convergence in the MSY calculations.}
  \item{endStrat}{\emph{numeric}: maximum fishing mortality for the MSY analysis.}
  \item{stepStrat}{\emph{numeric}: fishing mortality step size for the MSY analysis.}
  \item{delim}{\emph{character} that delimits the components of the \code{filename} prefix.}
  \item{clean}{\emph{logical}: if \code{TRUE}, clean all Awatea files before running an MPD analysis.}
  \item{locode}{\emph{logical}: if \code{TRUE}, source local code located at \code{codePath}.}
  \item{awateaPath}{\emph{character}: system path to the \code{Awatea.exe} executable (using path delimiter `/`).}
  \item{codePath}{\emph{character}: system path to \pkg{PBSawatea} source code (using path delimiter `/`).}
  \item{\dots}{additional arguments (not currently used for any purpose).}
}
\details{
  This function is primarily used to automate MPD reweightings and to perform 
  the MSY calculations.  MCMCs are better run from a command line console on 
  a supercomputer as they generally require \eqn{>12}h to complete.

  Once an MCMC has been created, the user can run various projections separately using \cr 
  \code{awatea -ind filename.ext -mceval} \cr on the command line.
}
\note{
  Be careful when setting the \code{clean} argument to \code{TRUE}. Some of the patterns for 
  matching to Awatea files might inadvertently clean non-related files. The patterns are currently set to: \cr
  \code{"^Awatea"}, \code{"^admodel"}, \code{"\\\\.pst$"}, \code{"\\\\.out$"}, \code{"\\\\.rpt$"}, \code{"\\\\.tmp$"}, \cr
  \code{"^variance$"}, \code{"^results.dat$"}, \code{"^likelihood.dat$"}
}
\section{Text Output Files Produced}{
  \code{Results.dat} contains an exhaustive listing of the maximum likelihood estimates
  for the model and derived parameters, a re-listing of some of the data and fixed
  parameters, and a listing of most of the prediction made. This includes numbers at age,
  fecundity, vulnerable biomass, survey trajectories, and so on. It is always placed in the
  same directory as the text input file. If there is an existing \code{Results.dat} file in that
  directory, it will be overwritten.

  Other output files generated during the parameter estimation process include \cr
  \code{Awatea.par} (contains the maximum likelihood estimates of the free parameters), \cr
  \code{Awatea.cor} (shows standard deviation and correlations between the estimated parameters), and \cr
  \code{Awatea.std} (standard deviation of the estimated parameters).

  Be aware that using argument \code{ADargs=list("-nohess")}
  will shorten the time to perform an MPD (useful if there will be multiple re-weightings), 
  but the output files \code{Awatea.cor} and \code{Awatea.std} will \bold{not} be produced.
}
\references{
Edwards, A.M., Haigh, R., and Starr, P.J. (2012a) 
  Stock assessment and recovery potential assessment for Yellowmouth Rockfish (\emph{Sebastes reedi}) along the Pacific coast of Canada.
  \emph{Canadian Science Advisory Secretariat, Research Document} \bold{2012/095}: iv + 188p.

Edwards, A. M., Haigh, R., and Starr, P.J. (2013a) 
  Pacific Ocean Perch (\emph{Sebastes alutus}) stock assessment for the north and west coasts of Haida Gwaii, British Columbia. 
  \emph{Fisheries and Oceans Canada. Canadian Science Advisory Secretariat Research Document} \bold{2013/092}: vi + 126p.

Edwards, A. M., Haigh, R., and P. J. Starr, P.J. (2013b) 
  Pacific Ocean Perch (\emph{Sebastes alutus}) stock assessment for the west coast of Vancouver Island, British Columbia. 
  \emph{Fisheries and Oceans Canada. Canadian Science Advisory Secretariat Research Document} \bold{2013/093}: vi + 135p.

Edwards, A.M., Starr, P.J., and Haigh, R. (2012b) 
  Stock assessment for Pacific Ocean Perch (\emph{Sebastes alutus}) in Queen Charlotte Sound, British Columbia.
  \emph{Canadian Science Advisory Secretariat, Research Document} \bold{2011/111}: viii + 172p.

Francis, R.I.C.C. (2011) 
  Data weighting in statistical fisheries stock assessment models. 
  \emph{Canadian Journal of Fisheries and Aquatic Sciences} \bold{68}: 1124-1138.

Hilborn, R., Maunder, M., Parma, A., Ernst, B. Payne, J., and Starr, P. (2003) 
  Coleraine: a generalized age-structured stock assessment model. 
  School of Aquatic and Fishery Sciences, University of Washington, 54p.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC.
}
\seealso{
  \code{\link{importRes}}, \code{\link{readAD}}, \code{\link{reweight}}, \code{\link{runSweave}}
}
\keyword{models}
\keyword{iteration}


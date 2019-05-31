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
  doMPD=FALSE,
  N.reweight = 1, A.reweight = 1, C.reweight = 1, cvpro=0,
  doMCMC=FALSE, 
  mcmc=1e6, mcsave=1e3, ADargs=NULL, verbose=FALSE, 
  doMSY=FALSE, 
  msyMaxIter=15000., msyTolConv=0.01, endStrat=0.301, stepStrat=0.001,
  delim="-", clean=FALSE, locode=FALSE, 
  awateaPath="C:/Users/haighr/Files/Projects/ADMB/Coleraine",
  codePath="C:/Users/haighr/Files/Projects/R/Develop/PBSawatea/Authors/Rcode/develop",
  ...
)
}
\arguments{
  \item{filename.ext}{\code{character} -- file name including its extension.}
  \item{wd}{\code{character} -- string specifying the working directory, if different from the current working directory.}
  \item{strSpp}{\code{character} -- string 3-letter code of the species.}
  \item{runNo}{\code{numeric} -- the model run number.}
  \item{rwtNo}{\code{numeric} -- the reweight number.}
  \item{doMPD}{\code{logical} -- if \code{TRUE}, perform an MPD analysis.}
  \item{N.reweight}{\code{numeric} -- the number of reweights to perform in the MPD analysis.}
  \item{A.reweight}{\code{numeric} -- abundance reweight method: 0=no reweight, 1=add pocess error, 2=SDNR.}
  \item{C.reweight}{\code{numeric} -- composition reweight method: 0=no reweight, 1=Francis (2011) mean age, 2=SDNR.}
  \item{cvpro}{\code{numeric} -- CV process error added to CV observation error: \eqn{c_t = \sqrt{c_o^2 + c_p^2}}; \cr
    if \code{cvpro} vector is shorter than \code{Nsurvey + Ncpue}, 
    it is repeated to match the number of survey series and CPUE series used to index abundance.}
  \item{doMCMC}{\code{logical} -- if \code{TRUE}, perform an MCMC analysis.}
  \item{mcmc}{\code{numeric} -- number of MCMC iterations to perform.}
  \item{mcsave}{\code{numeric} -- frequency of MCMC iterations to save.}
  \item{ADargs}{\code{character} -- additional arguments for a call to Awatea.}
  \item{verbose}{\code{logical} -- if \code{TRUE}, spew Awatea messages to the R console.}
  \item{doMSY}{\code{logical} -- if \code{TRUE}, perform an MSY analysis.}
  \item{msyMaxIter}{\code{numeric} -- maximum iterations for the MSY calculations.}
  \item{msyTolConv}{\code{numeric} -- tolerance for convergence in the MSY calculations.}
  \item{endStrat}{\code{numeric} -- maximum fishing mortality for the MSY analysis.}
  \item{stepStrat}{\code{numeric} -- fishing mortality step size for the MSY analysis.}
  \item{delim}{\code{character} that delimits the components of the \code{filename} prefix.}
  \item{clean}{\code{logical} -- if \code{TRUE}, clean all Awatea files before running an MPD analysis.}
  \item{locode}{\code{logical} -- if \code{TRUE}, source local code located at \code{codePath}.}
  \item{awateaPath}{\code{character} -- system path to the \code{Awatea.exe} executable (using path delimiter `/`).}
  \item{codePath}{\code{character} -- system path to \pkg{PBSawatea} source code (using path delimiter `/`).}
  \item{\dots}{additional arguments (not currently used for any purpose).}
}
\details{
  This function is primarily used to automate MPD reweightings and to perform 
  the MSY calculations.

  Reweighting of abundance data and composition data can now be performed independently of
  each other using the vectors \code{A.reweight} and \code{C.reweight}, respectively.
  This flexibility also means that a user could choose different reweight techniques
  for every reweight, but mixing techniques is not generally done.

  Typically, a user might perform one reweight and adjust both abundance and composition
  using added CV process error and Francis (2011) mean age adjustment of effective sample size:\cr
  \code{N.reweight=1, A.reweight=1, C.reweight=1, cvpro=0.2}

  Alternatively, a user might choose to reweight the composition data first and then the abundance data:\cr
  \code{N.reweight=2, A.reweight=c(0,1), C.reweight=c(1,0), cvpro=0.2}

  MCMCs are better run from a command line console on 
  a stand-alone computer as they generally require \eqn{>12}h to complete.
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
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2019-04-15}
}
\seealso{
  \code{\link{importRes}}, \code{\link{readAD}}, \code{\link{reweight}}, \code{\link{runSweave}}
}
\keyword{models}
\keyword{iteration}


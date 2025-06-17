\name{importProjRec}
\alias{importProjRec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Import Projected Recruitment
}
\description{
  Import the projected recruitments (actually, the values 
  are random normals \eqn{N(0,1)}).
}
\usage{
importProjRec(dir, info="", coda=FALSE, ngear=1, sigmaR, quiet=TRUE)
}
\arguments{
  \item{dir}{\code{character} -- directory where MCMC projections reside.}
  \item{info}{\code{character} -- user-supplied information, if desired.}
  \item{coda}{\code{logical} -- if \code{TRUE}, use the function \code{mcmc}
    in the package \pkg{coda} to generate MCMCs.}
  \item{ngear}{\code{numeric} -- number of commercial gears}
  \item{sigmaR}{\code{numeric} -- standard deviation of recruitment residuals.}
  \item{quiet}{\code{logical} -- if \code{TRUE}, print progress messages to the R console.}
}
\details{
  The values saved by the Awatea code are random normals \eqn{N(0,1)}, 
  which for a particular MCMC sample are the same for all the catch strategies.
}
\value{
  A list object comprising:
    \item{B}{data frame of spawning biomass (dim = MCMC samples by projected years)}
    \item{VB}{data frame of vulnerable biomass (dim = MCMC samples by projected years)}
    \item{Y}{data frame of yield (dim = MCMC samples by projected years)}
    \item{eps}{data frame of \eqn{\epsilon_t}{\epsilon[t]} (dim = MCMC samples by projected years)}
}
\author{
  \href{mailto:andrew.edwards@dfo-mpo.gc.ca}{Andrew Edwards}, Research Scientist, PBS, DFO, Nanaimo BC\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\note{
  The function \code{importProj} does not import recruitment residuals.

  This function grabs the \code{tempdev} values from Awatea, 
  which are just \eqn{N(0,1)} values, then multiplies them by \eqn{\sigma_R}{\sigma[R]}
  to yield \eqn{\epsilon_t \sim N(0,\sigma_R^2)}{\epsilon[t] ~ N(0,\sigma[R]^2)}.

  The parameter value for \eqn{\sigma_R}{sigma[R]} can be found in \cr
  \code{currentRes$extra$residuals$p_log_RecDev[6]}.
}
\seealso{
  \code{\link{importProj}} \cr
  \pkg{coda}: \code{\link[coda]{mcmc}}
}
\keyword{ts}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\concept{Projection}

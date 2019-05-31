\name{splineCPUE}
\alias{splineCPUE}
\title{
  Spline-smooth CPUE Indices
}
\description{
  Fit spline curves through CPUE data to determine optimal degrees of freedom (DF)
  -- balance between rigorously fitting indices while not removing the majority
  of the signal -- and calculate CV process error from residual sum of squares
  (RSS) at optimal DF.
}
\usage{
splineCPUE(dat, ndf=50, png=FALSE, pngres=400, PIN=c(8,8))
}
\arguments{
  \item{dat}{\code{numeric} -- CPUE index series as a vector of length \eqn{N}.}
  \item{ndf}{\code{numeric} -- number of degrees of freedom between 2 and \eqn{N}.}
  \item{png}{\code{logical} -- if \code{TRUE}, send the four-panel figure to a \code{.png} file.}
  \item{pngres}{\code{numeric} -- resolution of \code{.png} file in pixels per inch.}
  \item{PIN}{\code{numeric} -- size (width, height) of the \code{.png} file.}
}
\details{
  The code uses a procedure for estimating process error \eqn{c_\mathrm{p}} to
  add to commercial CPUE using a spline-smoother analysis. 
  Francis (2011, citing Clark and Hare 2006) recommends using a smoothing
  function to determine the appropriate level of process error to add to CPUE
  data, with the goal of finding a balance between rigorously fitting the 
  indices while not removing the majority of the signal in the data.
  An arbitrary sequence of length \code{'ndf'}, comprising degrees of freedom
  (DF, \eqn{\nu_i}), where \eqn{i=2,...,N} and \eqn{N}= number of CPUE values \eqn{U_t} from
  \eqn{t=\mathrm{yr}_1,...,\mathrm{yr}_N}, is used to fit the CPUE data with a spline smoother.
  At \eqn{i=N}, the spline curve fits the data perfectly and the residual sum of 
  squares (RSS, \eqn{\rho_N}) is 0.
  Using spline fits across a range of trial DF \eqn{\nu_i}, values of RSS \eqn{\rho_i}
  form a logistic-type curve with an inflection point at \eqn{i=k}.
  The difference between point estimates of \eqn{\rho_i} (proxy for the slope 
  \eqn{\delta_i}) yields a concave curve with a minimum \eqn{\delta_i}, which 
  occurs close to the inflection point \eqn{k}.
  At the inflection point \eqn{k}, \eqn{\nu_k} corresponds to \eqn{\rho_k}, which is 
  converted to \eqn{c_\mathrm{p}} using:

  \deqn{c_\mathrm{p} = \sqrt{\frac{\rho_k}{N-2}}~~~{\left[ \frac{1}{N} \sum\limits_{t=\mathrm{yr}_1}^{\mathrm{yr}_N} U_t \right]}^{-1}}{%
    c[p] = \sqrt{\rho[k]/(N-2)} / ( (1/N) \Sigma_{t=yr[1]}^{yr[N]} U[t] }
}
\value{
  Returns a list containg the optimal degrees of freedom (DF), the residual sum
  of squares (RSS) at the lowest first derivative of RSS, and the calculated 
  process error at this RSS (see \bold{Details}).
}
\references{
  Clark, W.G. and Hare, S.R. 2006.
  \href{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.170.6113&rep=rep1&type=pdf}{Assessment and management of Pacific halibut: data, methods, and policy.}
  Sci. Rep. 83, International Pacific Halibut Commission, Seattle, WA.

  Francis, R.I.C.C. 2011.
  \href{http://dx.doi.org/10.1139/f2011-025}{Data weighting in statistical fisheries stock assessment models}.
  Can. J. Fish. Aquat. Sci. 68(6): 1124-1138.

  Starr, P.J. and Haigh, R. 2019, in review.
  Widow Rockfish (\emph{Sebastes entomelas}) stock assessment for British Columbia in 2019. 
  DFO Can. Sci. Advis. Sec. Res. Doc. 2019/nnn. vi + 219 p.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2019-05-28}
}
\note{
  This procedure was first used in the 2019 BC Widow Rockfish stock assessment (Starr and Haigh, in review).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
\keyword{optimize}% __ONLY ONE__ keyword per line


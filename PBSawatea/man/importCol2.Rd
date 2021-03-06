\name{importCol2}
\alias{importCol2}
\title{
  Import Coleraine Model Results (AME version)
}
\description{
  Import Coleraine model results from \code{.res} file, and 
  rearrange into a standard format suitable for plotting. 
}
\usage{
importCol2(res.file, info="", Dev=FALSE, CPUE=FALSE, Survey=FALSE, 
    CAc=FALSE, CAs=FALSE, CLc=FALSE, CLs=FALSE, LA=FALSE, 
    quiet=TRUE, extra=TRUE)
}
\arguments{
  \item{res.file}{name of Coleraine model results file to import.}
  \item{info}{optional string containing information to store with model results.}
  \item{Dev}{logical: whether recruitment deviates were estimated in model.}
  \item{CPUE}{logical: whether model was fitted to catch-per-unit-effort data.}
  \item{Survey}{logical: whether model was fitted to survey abundance index data.}
  \item{CAc}{logical: whether model was fitted to commercial catch-at-age data.}
  \item{CAs}{logical: whether model was fitted to survey catch-at-age data.}
  \item{CLc}{logical: whether model was fitted to commercial catch-at-length data.}
  \item{CLs}{logical: whether model was fitted to survey catch-at-length data.}
  \item{LA}{logical: whether model was fitted to length-at-age data.}
  \item{quiet}{logical: whether to report progress while parsing file.}
  \item{extra}{logical: if \code{TRUE}, import likelihoods, parameters, priors, and recruitment residuals.}
}
\details{
  This function was modified from the original \code{importCol} function
  in the \pkg{scape} package to grab extra data.
}
\value{
  A list of class \code{scape} containing at least \code{N}, \code{B},
  and \code{Sel}. The other elements may or may not be included in the
  list, depending on how \code{importRes} was called:
  \item{N}{predicted numbers at age}
  \item{B}{predicted biomass, recruitment, and observed landings (year
    things)}
  \item{Sel}{predicted selectivity and observed maturity (age things)}
  \item{Dev}{predicted recruitment deviates from the stock-recruitment
    curve}
  \item{CPUE, Survey}{commercial and survey abundance index and fit}
  \item{CAc, CAs}{commercial and survey C@A (catch at age) and fit}
  \item{CLc, CLs}{commercial and survey C@L (catch at length) and fit}
  \item{LA}{observed L@A and fit}
}
\references{
  Hilborn, R., Maunder, M., Parma, A., Ernst, B., Payne, J., and Starr, P. (2003) 
  Coleraine: A generalized age-structured stock assessment model. User's manual version 2.0. 
  \emph{University of Washington Report} \bold{SAFS-UW-0116}. Available at: \cr
  \url{http://fish.washington.edu/research/coleraine/pdf/coleraine.pdf}.
}
\note{
  This import function is implemented for the Coleraine statistical
  catch-at-age software, and can serve as a template for \pkg{scape}
  users who would like to implement import functions for specific stock
  assessment software.

  The functions \code{ll} (package \pkg{gdata}) and \code{head} are
  recommended for browsing model results, \cr 
  \emph{e.g.} \code{ll(x.cod); ll(x.cod$N); head(x.cod$N)}.
}
\seealso{
  \code{\link{importRes}}, \code{\link[utils]{read.table}}, \code{\link[base]{readLines}}, 
  and \code{\link[base]{scan}} to import any data. \cr
  \code{\link[scape]{scape-package}} gives an overview of the package \pkg{scape}.
}
\keyword{file}
\keyword{interface}


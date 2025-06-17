\name{importRes}
\alias{importRes}
\title{
  Import Coleraine Model Results (RH version)
}
\description{
  Import Coleraine model results from \code{.res} file, and 
  rearrange into a standard format suitable for plotting. 
}
\usage{
importRes(res.file, info="", Dev=FALSE, CPUE=FALSE, Survey=FALSE,
    CAc=FALSE, CAs=FALSE, CLc=FALSE, CLs=FALSE, LA=FALSE, quiet=TRUE,
    extra=TRUE, sep=" ")
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
  \item{sep}{the field separator character (usually \code{" "} or \code{"\t"}). }
}
\details{
  This function was modified from the original \code{importCol} function
  in the \pkg{scape} package to grab extra data and to deal with anomalous 
  characters in Coleraine results files.
}
\value{
  A list of class \code{list} containing at least \code{N}, \code{B},
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
  \href{https://www.researchgate.net/publication/33515338_COLERAINE_A_Generalized_Age-Structured_Stock_Assessment_Model-User\%27s_Manual_Version_20}{%
  Coleraine: A generalized age-structured stock assessment model. User's manual version 2.0}. 
  \emph{University of Washington Report} \bold{SAFS-UW-0116}.
}
\author{
Original author of \code{scape::importCol} :\cr
  \href{mailto:arnim@spc.int}{Arni Magnusson}, Senior Fisheries Scientist at SPC\cr
  Intergovernmental organisation for sustainable development of the Pacific islands, Noum\ifelse{latex}{\out{\'{e}}}{\ifelse{html}{\out{&eacute;}}{e}}a , New Caledonia
%%https://stackoverflow.com/questions/71021330/how-to-insert-accented-characters-in-ascii-encoded-rd-files

Modified version \code{PBSawatea::importRes} :\cr
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
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
  \code{\link[PBSawatea]{importCor}},
  \code{\link[PBSawatea]{importPar}},
  \code{\link[PBSawatea]{importStd}},
  \code{\link[PBSawatea]{runADMB}},
  \code{\link[PBSawatea]{readAD}},
  \code{\link[PBSawatea]{reweight}}\cr
  Import data:\cr
  \code{\link[utils]{read.table}},
  \code{\link[base]{readLines}},
  \code{\link[base]{scan}}\cr
  Overview of \pkg{scape} package:\cr
  \code{\link[scape]{scape-package}}
}
\keyword{file}
\keyword{interface}


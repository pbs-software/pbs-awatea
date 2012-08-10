	\name{importLik}
\alias{importLik}
\title{
  Import Awatea Likelihoods
}
\description{
  Import an Awatea file of the likelihoods from an MPD estimation.
}
\usage{
importLik(lik.file)
}
\arguments{
  \item{lik.file}{likelihoods file (usually \code{likelihood.dat})}
}
\details{
  Extracts the data from character vectors and makes various data objects.
}
\value{
  The output is a list object with the following components:
  \item{lik}{character vector representing the file line-by-line}
  \item{Total_likelihood}{values of the final likelihood (objective function value)}
  \item{Worst_gradient}{worst gradient?}
  \item{CPUE}{commercial CPUE likelihood(s)}
  \item{Survey_Index}{survey index likelihood(s)}
  \item{CA_Commercial}{commercial catch-at-age likelihood(s)}
  \item{CL_Commercial}{commercial catch-at-length likelihood(s)}
  \item{CA_survey}{survey catch-at-age likelihood(s)}
  \item{CL_no_sex_data_survey}{likelihood(s) for catch-at-length from surveys with no sex data?}
  \item{CL_data_survey}{likelihood(s) for catch-at-length from surveys with sex data?}
  \item{Von_B_Likelihood}{von Bertalanffy likelihood}
  \item{Priors}{priors?}
  \item{Penalty_for_U}{penalty for U?}
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link[PBSawatea]{importCor}}, \code{\link[PBSawatea]{importPar}},
  \code{\link[PBSawatea]{importRes}}, \code{\link[PBSawatea]{importStd}}
}
\keyword{file}
\keyword{interface}


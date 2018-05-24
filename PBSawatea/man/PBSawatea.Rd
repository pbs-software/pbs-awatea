\name{PBSawatea}
\docType{package}
\alias{PBSawatea}
\alias{PBSawatea-package}
\title{PBS Awatea}

\description{
  \pkg{PBSawatea} contains the code used for the modelling of populations of 
  Pacific Ocean Perach (\emph{Sebastes alutus}) and Yellowmouth Rockfish 
  (\emph{S. reedi}) along the British Columbia (BC) coast.

  Implementation is done using a modified version of the Coleraine statistical catch-at-age
  software (Hilborn \emph{et al.} 2003) called Awatea (Alan Hicks, NOAA, pers. comm.). Awatea
  is a platform for implementing the AD (Automatic Differentiation) Model Builder software
  (Otter Research 1999), which provides (a) maximum posterior density estimates using a
  function minimiser and automatic differentiation, and (b) an approximation of the posterior
  distribution of the parameters using the Markov Chain Monte Carlo (MCMC) method,
  specifically using the Hastings-Metropolis algorithm (Gelman \emph{et al.} 2004).

  Running of Awatea is streamlined using code written in R (R Development Core Team 2009),
  rather than the original Microsoft Excel implementation. Figures and tables of output are
  automatically produced through R using code adapted from the R packages \pkg{scape} 
  (Magnusson 2009) and \pkg{plotMCMC} (Magnusson and Stewart 2007). We use the R function 
  \code{Sweave} (Leisch 2008) in the package \pkg{utils} to automatically collate, 
  via LATEX, the large amount of figures and tables into a single portable document file
  (\code{.pdf}) for each model run.
  
  We  provide master \code{Sweave} files used in folder \code{../library/PBSawatea/snw}
  to build the \code{.pdf} document. The user must copy these to a local working directory 
  if they are not already there.\cr
}
\references{
  Gelman, A., Carlin, J.B., Stern, H.S, and Rubin, D.B. (2004)
  Bayesian data analysis, 2nd edition. Chapman and Hall/CRC, New York, 668 p.

  Hilborn, R., Maunder, M., Parma, A., Ernst, B. Payne, J., and Starr, P. (2003)
  Coleraine: a generalized age-structured stock assessment model. 
  \emph{School of Aquatic and Fishery Sciences}, University of Washington, 54 p.

  Leisch, F. (2008) Sweave, R package.

  Magnusson, A. (2009)
  Scape -- statistical catch-at-age plotting environment, R package.

  Magnusson, A. and Stewart, I. (2007)
  MCMCscape -- MCMC diagnostic plots. R package.

  Otter Research Ltd. (1999)
  An introduction to AD Model Builder for use nonlinear modeling and statistics.
  Otter Research Ltd., British Columbia. 194 p.

  R Development Core Team (2011)
  R: A language and environment for statistical computing. 
  R Foundation for Statistical Computing, Vienna, Austria. \cr
  ISBN 3-900051-07-0 \cr
  \url{https://www.r-project.org/}
}
\keyword{package}

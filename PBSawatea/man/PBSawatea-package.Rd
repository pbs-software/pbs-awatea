\name{PBSawatea-package}
\docType{package}
\alias{PBSawatea-package}
\alias{PBSawatea}

\title{PBSawatea -- Tools for Running Awatea and Visualizing the Results}

\description{
  \pkg{PBSawatea} contains R code to control the input to and output from the 
  modelling of rockfish populations (e.g., Pacific Ocean Perch, \emph{Sebastes alutus})
  that occur along the British Columbia (BC) coast.

  The model implementation uses a modified version of the Coleraine statistical catch-at-age
  software (Hilborn \emph{et al.} 2003) called Awatea (Allan Hicks, IPHC, pers. comm.).
  Awatea template code is compiled by Automatic Differentiation Model Builder
  (ADMB) software (Otter Research 1999) to a 64-bit binary called \sQuote{\code{awatea.exe}}.
  Running the binary file, which contains ADMB libraries, provides
  (a) maximum posterior density estimates using a function minimiser and automatic differentiation, and 
  (b) an approximation of the posterior distribution of the parameters using the Markov Chain Monte Carlo (MCMC)
  method, specifically using the Hastings-Metropolis algorithm (Gelman et al. 2004).

  The running of Awatea is streamlined, using code written in R (R Development Core Team 2018)
  rather than manually manipulating the original Microsoft Excel implementation.
  Figures and tables of output are automatically produced through R using code adapted
  from the R packages \pkg{scape} (Magnusson 2009) and \pkg{plotMCMC} (Magnusson & Stewart 2007).
  We use the R function \code{Sweave} (Leisch 2008) in the package \pkg{utils} to 
  automatically collate, via LaTeX (specifically pdflatex), the large amount of figures
  and tables into a single portable document file (\code{.pdf}) for each model run.
  
  We  provide master \code{Sweave} files in the folder \code{../library/PBSawatea/snw}
  to build the \code{.pdf} document. The code copies these to a local working directory 
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

  R Development Core Team (2018)
  R: A language and environment for statistical computing. 
  href{https://www.r-project.org/}{R Foundation for Statistical Computing}, Vienna, Austria.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-25}
}
\keyword{package}


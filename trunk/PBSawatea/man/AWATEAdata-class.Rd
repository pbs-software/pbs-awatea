\name{AWATEAdata-class}
\docType{class}
\alias{AWATEAdata-class}
\alias{fix,AWATEAdata-method}
\alias{reweight,AWATEAdata-method}
\alias{view,AWATEAdata-method}
\alias{write,AWATEAdata-method}

\title{S4: AWATEA Data Class}
\description{
	The set of functions described here provide an easy method of dealing with a 
	complicated scheme for manipulating data when using the catch-at-age model 
	called \pkg{Awatea}, a variant of \pkg{Coleraine}. 
	\code{Awatea} is maintained by Allan Hicks (University of Washington) and Paul Starr
	(Canadian Groundfish Research and Conservation Society).
	
	To use \pkg{PBSawatea}'s management functions, the function \code{runADMB} calls \code{readAD} 
	which creates an \code{AWATEAdata} object. The \code{AWATEAdata} object contains a distinct 
	R environment where data are stored and accessed by various package functions.

}
\section{Objects from the Class}{
  This S4 object can be created by calls of the form: 

  \code{dat = new("AWATEAdata", txtnam=character(), input=character(), vlst=list(),} \cr
  \code{      dnam=character(), nvars=numeric(), vdesc=character(), vars=list(), } \cr
  \code{      gcomm=character(), vcomm=character(), resdat=list() , likdat=list(), } \cr
  \code{      pardat=list(), stddat=list(), cordat=list(), evadat=list(), reweight=list()) } \cr

  The function \code{readAD} populates the first nine slots of the S4 object by reading in an \code{Awatea} data
  file from the user's current working directory. If \code{Awatea} output files exist (with the same prefix 
  as the input file and suffixes = \{\code{.res}, \code{.lik}, \code{.par}, \code{.std}, \code{.cor}, \code{.eva}\}),
  then slots \code{resdat}, \code{likdat}, \code{pardat}, \code{stddat}, \code{cordat}, and \code{evadat}, respectively,
  will also be populated.
}
\section{Slots}{
  \describe{
    \item{\code{txtnam}}{the file name (including extension) of the \code{Awatea} input data file. }
    \item{\code{input}}{a string vector comprising the line-by-line information from the data file. }
    \item{\code{vlst}}{a list named by line number (e.g. \code{"L001:..."}) of comments and data. }
    \item{\code{dnam}}{a vector of data descriptors referenced by numeric line number. }
    \item{\code{nvars}}{a numeric scalar indicating the number of input variables. }
    \item{\code{vdesc}}{as for \code{dnam} but indexed by variable identifier (e.g. \code{"v001"}).}
    \item{\code{vars}}{variable values indexed by variable identifier. }
    \item{\code{gcomm}}{commented lines indexed by numeric line number. }
    \item{\code{vcomm}}{value descriptors indexed by numeric line number.}
    \item{\code{resdat}}{list of Awatea results from the file \code{results.dat} after calling \code{importRes}.}
    \item{\code{likdat}}{list of MPD likelihoods from the file \code{likelihood.dat} after calling \code{importLik}.}
    \item{\code{pardat}}{list of parameter values from the file \code{Awatea.par} after calling \code{importPar}.}
    \item{\code{stddat}}{list of estimated parameter values from the file \code{Awatea.std} after calling \code{importStd}.}
    \item{\code{cordat}}{list of correlation objects from the file \code{Awatea.cor} after calling \code{importCor}.}
    \item{\code{evadat}}{eigenvalues of the Hessian from the file \code{Awatea.eva} after calling \code{importEva}.}
    \item{\code{reweight}}{list of re-weighted abundance and composition data (see \code{\link[PBSawatea]{reweight}} for details).}
  }
}
\section{Methods}{
  \describe{
    \item{\code{fix, signature(x="AWATEAdata")}}{:  replace data elements of an input list }
    \item{\code{reweight, signature(x="AWATEAdata")}}{:  reweight abundance and composition data }
    \item{\code{view, signature(x="AWATEAdata")}}{:  view the basic input list }
    \item{\code{write, signature(x="AWATEAdata")}}{:  write a new input data file }
  }
}
\details{
  An \code{Awatea} data file can be loaded into an \code{AWATEAdata} object using \code{readAD},
  and if a complementary results file is available, the user can run \code{reweight} for a single
  re-weighting of the abundance data and the composition data.
  
  Alternatively, starting with a single input data file, a user can run the function \code{runADMB} 
  to perform multiple (iterative) re-weightings automatically.
}
\note{
  Some of the output values in \code{results.dat}, \code{likelihood.dat}, \code{Awatea.par},
  \code{Awatea.std}, and \code{Awatea.cor} contain redundant information.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}

\seealso{
  \code{\link[PBSawatea]{readAD}} for loading the ADMB data file

  \code{\link[PBSawatea]{importRes}} for reading in the file \code{results.dat}

  \code{\link[PBSawatea]{importLik}} for reading in the file \code{likelihood.dat}

  \code{\link[PBSawatea]{importPar}} for reading in the file \code{Awatea.par}

  \code{\link[PBSawatea]{importStd}} for reading in the file \code{Awatea.std}

  \code{\link[PBSawatea]{importCor}} for reading in the file \code{Awatea.cor}

  \code{\link[PBSawatea]{importEva}} for reading in the file \code{Awatea.eva}

  \code{\link[PBSawatea]{reweight}} for re-weighting abundance and composition data

  \code{\link[PBSawatea]{runADMB}} for running the \code{Awatea} model through multiple iterations.
}
\keyword{classes}
\keyword{methods}


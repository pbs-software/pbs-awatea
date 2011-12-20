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
	\code{Awatea} is maintained by Alan Hicks (University of Washington) and Paul Starr
	(Canadian Groundfish Research and Conservation Society).
	
	To use \pkg{PBSawatea}'s management functions, the function \code{runADMB} calls \code{readAD} 
	which creates an \code{AWATEAdata} object. The \code{AWATEAdata} object contains a distinct 
	R environment where data are stored and accessed by various package functions.

}
\section{Objects from the Class}{
  This S4 object can be created by calls of the form: 

  \code{dat = new("AWATEAdata", txtnam=character(), input=character(), vlst=list(),} \cr
  \code{      dnam=character(), nvars=numeric(), vdesc=character(), vars=list(), } \cr
  \code{      gcomm=character(), vcomm=character(), output=list() , reweight=list()) }


  The function \code{readAD} populates the first nine slots of the S4 object by reading in an \code{Awatea} data
  file from the user's current working directory. If an \code{Awatea} results file (with the same prefix 
  as the input file and the suffix = \code{.res}), then slot 10 \code{output} will also be populated.
}
\section{Slots}{
  \describe{
    \item{\code{txtnam}}{: the file name (including extension) of the \code{Awatea} input data file. }
    \item{\code{input}}{: a string vector comprising the line-by-line information from the data file. }
    \item{\code{vlst}}{: a list named by line number (e.g. \code{"L001:..."}) of comments and data. }
    \item{\code{dnam}}{: a vector of data descriptors referenced by numeric line number. }
    \item{\code{nvars}}{: a numeric scalar indicating the number of input variables. }
    \item{\code{vdesc}}{: as for \code{dnam} but indexed by variable identifier (e.g. \code{"v001"}).}
    \item{\code{vars}}{: variable values indexed by variable identifier. }
    \item{\code{gcomm}}{: commented lines indexed by numeric line number. }
    \item{\code{vcomm}}{: value descriptors indexed by numeric line number.}
    \item{\code{output}}{: list of variable values from the \code{Awatea} results file using \code{importRes}.}
    \item{\code{reweight}}{: list of re-weighted abundance and composition data (see \code{\link[PBSawatea]{reweight}} for details).}
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
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}

\seealso{
  \code{\link[PBSawatea]{readAD}} for loading the ADMB data file

  \code{\link[PBSawatea]{importRes}} for reading in an \code{Awatea} results file

  \code{\link[PBSawatea]{reweight}} for re-weighting abundance and composition data

  \code{\link[PBSawatea]{runADMB}} for running the \code{Awatea} model through multiple iterations.
}
\keyword{classes}
\keyword{methods}


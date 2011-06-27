\name{readAD}
\alias{readAD}
\title{
  Read ADMB Input and Create AWATEA Class Object
}
\description{
  Read the ADMB input file and create an \code{AWATEA} class object.
}
\usage{
readAD(txt)
}
\arguments{
  \item{txt}{string name of an Awatea input file.}
}
\details{
  The Awatea input file contains headers (lines prefixed with \code{"#"}) and data
  that are read sequentially into the model by the binary executable \code{Awatea.exe}.
}
\value{
  An \code{AWATEA} class cobject with the slots: \cr
  \code{txtnam.....}{character: name of the input file,} \cr
  \code{input......}{character: vector of strings that are the lines of the input file,} \cr
  \code{vlst.......}{list: each line of the input file with a label specifying line number and} \cr
  \code{...........}{indicating whether the line is a Comment or Data,} \cr
  \code{dnam.......}{character: vector of strings specifying data contents labelled by line number,} \cr
  \code{nvars......}{numeric: number of data variables,} \cr
  \code{vdesc......}{character: vector of strings specifying data contents labelled by} \cr
  \code{...........}{variable number (\emph{e.g.}, \code{v001}),} \cr
  \code{vars.......}{list: numeric values of data labelled by variable number,} \cr
  \code{gcomm......}{character: vector of comments labelled by line number,} \cr
  \code{vcomm......}{character: vector of variable names labelled by line number,} \cr
  \code{output.....}{list: Awatea results file imported by function \code{\link{importRes}},} \cr
  \code{reweight...}{list: empty (later populated by function \code{\link{reweight}}).}
}
\references{
  Hilborn, R., Maunder, M., Parma, A., Ernst, B. Payne, J., and Starr, P. (2003)
  Coleraine: a generalized age-structured stock assessment model. 
  School of Aquatic and Fishery Sciences, University of Washington, 54 p.
}
\author{
  Rowan Haigh, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC.
}
\seealso{
  \code{\link{runADMB}}, \code{\link{reweight}}
}
\keyword{models}
\keyword{methods}


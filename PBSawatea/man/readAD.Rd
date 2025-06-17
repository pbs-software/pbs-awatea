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
  \item{txt}{\code{character} -- name of an Awatea input file}
}
\details{
  The Awatea input file contains headers (lines prefixed with \code{"#"}) and data
  that are read sequentially into the model by the binary executable \code{Awatea.exe}.
}
\value{
  An \code{AWATEA} class cobject with the slots:
  \tabular{ll}{
    \code{txtnam}   \tab character: name of the input file\cr
    \code{input}    \tab character: vector of strings that are the lines of the input file\cr
    \code{vlst}     \tab list: each line of the input file with a label specifying line number and\cr
                    \tab indicating whether the line is a Comment or Data\cr
    \code{dnam}     \tab character: vector of strings specifying data contents labelled by line number\cr
    \code{nvars}    \tab numeric: number of data variables\cr
    \code{vdesc}    \tab character: vector of strings specifying data contents labelled by\cr
                    \tab variable number (\emph{e.g.}, \code{v001})\cr
    \code{vars}     \tab list: numeric values of data labelled by variable number\cr
    \code{gcomm}    \tab character: vector of comments labelled by line number\cr
    \code{vcomm}    \tab character: vector of variable names labelled by line number\cr
    \code{output}   \tab list: Awatea results file imported by function \code{\link{importRes}}\cr
    \code{reweight} \tab list: empty (later populated by function \code{\link{reweight}})
  }
}
\references{
  Hilborn, R., Maunder, M., Parma, A., Ernst, B. Payne, J., and Starr, P. (2003)
  Coleraine: a generalized age-structured stock assessment model. 
  School of Aquatic and Fishery Sciences, University of Washington, 54 p.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-17}
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{runADMB}},
  \code{\link[PBSawatea]{reweight}}
}
\keyword{models}
\keyword{methods}
\concept{runADMB.r}


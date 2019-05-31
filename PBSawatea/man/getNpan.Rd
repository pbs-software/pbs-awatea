\name{getNpan}
\alias{getNpan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get Number of Current Panel
}
\description{
  Get the number of the current panel when plotting a multi-panel plot.
  Useful when plots are called from an \code{'sapply'} or when functions
  call a panel subroutine.
}
\usage{
getNpan()
}
\details{
  Locates the panel based on \code{'par()$mgp'}. Basically, the function is a wrapper for:
  \preformatted{
  mfg=par()$mfg
  mfg[2]+(mfg[1]-1)*mfg[4]
  }
}
\value{
  The current panel, assuming plots are executed by row.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2019-05-29}
}
\note{
  May not work if you are using \code{'par(mfcol=c(r,c))'}.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

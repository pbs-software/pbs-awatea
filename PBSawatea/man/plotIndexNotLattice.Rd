\name{plotIndexNotLattice}
\alias{plotIndexNotLattice}
\title{
  Plot Survey Indices
}
\description{
  Plot survey index series with error bars in a variety of ways (e.g., grouped, individually).
}
\usage{
plotIndexNotLattice(obj, main="", save=NULL,
   bar=1.96, ssnames=paste("Ser",1:9,sep=""),
   ptypes=tcall(PBSawatea)$ptype, pngres=400, lang=c("e","f"), ...)
}
\arguments{
  \item{obj}{\code{list} -- the MPD object created by the \code{importRes} function (e.g., \code{currentRes}).}
  \item{main}{\code{character} -- title for figure}
  \item{save}{\bold{not used}}
  \item{bar}{\code{numeric} -- standard deviation of the normal distribution 
    (1.96 is the approximate value of the 97.5 percentile point).}
  \item{ssnames}{\code{character} -- survey series names for figure labelling and image creation.}
  \item{ptypes}{\code{character} -- plot types for the figures: \code{"eps"} = encapsulated postscript, \code{"png"} = portable network graphics.}
  \item{pngres}{\code{numeric} -- resolution (pixels/inch) for the \code{.png} file.}
  \item{lang}{\code{character} -- a vector of letters that denote the language for output:
    currently only \code{"e"} (english) and \code{"f"} (french).}
  \item{\dots}{\bold{not currently used}}
}
\value{
  Four series of encapsulated postscript (and portable network graphics) files: \cr
  \code{survIndSer.eps....}{each survey panel focuses on the years of the survey;} \cr
  \code{survIndSer2.eps...}{each panel uses a fixed set of years that span all surveys;} \cr
  \code{survIndSer3.eps...}{one panel showing all series (incl. CPUE if used) normalised to their means;} \cr
  \code{survIndSer4-n.eps.}{seperate files for each index series \emph{n}.}
}
\author{
  Andrew M. Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC

  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-07-26}
}
\seealso{
  \code{\link{plt.idx}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


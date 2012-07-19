\name{plotIndexNotLattice}
\alias{plotIndexNotLattice}
\title{
  Plot Survey Indices
}
\description{
  Plot index series with error bars. Create postscript files automatically.
}
\usage{
plotIndexNotLattice(obj, objCPUE, main="", save=NULL, bar=1.96,
    ssnames=paste("Ser",1:9, sep=""), ...)
}
\arguments{
  \item{obj}{data frame of survey indices from Awatea's results file
    (\emph{e.g.}, \code{currentRes$Survey}).}
  \item{objCPUE}{data frame of CPUE indices from Awatea's results file
    (\emph{e.g.}, \code{currentRes$CPUE}).}
  \item{main}{title for figure}
  \item{save}{\bold{not used}}
  \item{bar}{standard deviation of the normal distribution 
    (1.96 is the approximate value of the 97.5 percentile point).}
  \item{ssnames}{survey series names for figure labelling and image creation.}
  \item{\dots}{\bold{not used}}
}
\value{
  Four postscript files: \cr
  \code{survIndSer.eps....}{each survey panel focuses on the years of the survey;} \cr
  \code{survIndSer2.eps...}{each panel uses a fixed set of years that span all surveys;} \cr
  \code{survIndSer3.eps...}{one panel showing all series normalised to their means;} \cr
  \code{survIndSer4.eps...}{compares normalised GIG series with CPUE series.}
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link{plt.idx}}
}
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


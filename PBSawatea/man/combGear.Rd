\name{combGear}
\alias{combGear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Combine Gear Values
}
\description{
  Combine MCMC values from multiple gears by year -- useful when data input specifies multiple gear types (fisheries).
}
\usage{
combGear(dat, fn = function(x) {
    sum(x, na.rm = TRUE)
})
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{\code{matrix|data.frame} -- usually harvest rates or vulnerable biomass.}
  \item{fn}{\code{function} -- summary function to apply across rows and combine columns of the same year.}
}
\details{
  If multiple gear types (fisheries) are specified in the input data file, MCMC outputs
  of harvest rate and vulnerable biomass have column headers like: 
  \code{c("1940_1", "1940_2", ..., "2018_1", "2088_2")}.
  The corresponding projections, however, only have one value per year.
  
  Typically, a user might want to represent vulnerable biomass by adding values in the same year.
  Harvest rates are more logically averaged, preferably using the geometric mean.
  However, users can supply their own functions for summarising annual values, should they wish.
  Currently only used in compiling the \sQuote{Model Results} for offshore rockfish stock assessments.
}
\value{
  The input data frame but with columns reduced by 1/Ngear and annual values summarised by the user's function.
  If the input data only has one column per year, no summary action is applied.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: remote office, Vancouver BC\cr
  Last modified \code{Rd: 2025-06-16}
}
\seealso{
In package \pkg{PBSawatea}:\cr
  \code{\link[PBSawatea]{plotSnail}},
  \code{\link[PBSawatea]{plotTraj}}
}	
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
\keyword{arith}% __ONLY ONE__ keyword per line

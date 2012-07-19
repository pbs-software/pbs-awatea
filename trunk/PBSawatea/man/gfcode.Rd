\name{gfcode}
\docType{data}
\alias{gfcode}
\title{Code: Groundfish Species Codes and Names}
\description{
  Data frame of groundfish species described primarily by Hart codes. The information resides in 
  a data table called \code{C_Species} in the relational database \code{PacHarvest}.
}
\usage{data(gfcode)}
\format{
  Data frame comprising 549 rows (species) and 5 columns:
  \tabular{ll}{
    \code{code}    \tab Hart code \cr
    \code{latin}   \tab scientific (latin) name \cr
    \code{name}    \tab common species name \cr
    \code{code2}   \tab 2-letter code if available \cr
    \code{code3}   \tab 3-letter code if available
  }
}
\details{
  The data frame comprises 5,49 rows (species) and 5 columns. Hart codes for fish 
  species essentially reference page numbers in Hart (1973).

  Quick reference to species name:  \code{species["424",]$name} yields \code{Quillback rockfish}.

  Quick reference to latin name:  \code{species["424",]$latin} yields \code{Sebastes maliger}.
}
\source{
  Norm Olsen, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\references{
  Hart, J.L. (1973) Pacific Fishes of Canada. \emph{Bulletin of the Fisheries Research 
  Board of Canada} \bold{180}, 740 pp.
}
\keyword{datasets}

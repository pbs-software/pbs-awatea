\name{out.pmTables}
\alias{out.pmTables}
\title{
  Write Decision Tables to Comma-Delimited Files
}
\description{
  Write decision tables to comma-delimited text files (\code{.csv}).
}
\usage{
out.pmTables(obj, fileName="pm", dec=3)
}
\arguments{
  \item{obj}{list object containing tables (matrices or data frames).}
  \item{fileName}{prefix for output file names.}
  \item{dec}{number of decimal places to retain.}
}
\value{
  Comma-delimited text files (\code{.csv}).
}
\author{
  Andrew Edwards, Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC
}
\seealso{
  \code{\link[utils]{write.table}}, \code{\link[PBSmodelling]{writeList}}
}
\keyword{IO}
\keyword{file}


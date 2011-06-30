\name{get.resFile}
\alias{get.resFile}
\title{
  Get Awatea Results Files for Menu
}
\description{
  A function that retrieves the names of Awatea results files (\code{.res$}) for
  use in the \code{mainMenu} command. When choice is made, the function loads the 
  results file and assigns it to the global environment as \code{currentRes}.
}
\usage{
get.resFile(resFile=NULL)
}
\arguments{
  \item{resFile}{supposedly the name of a results file, but the code
    suggests that argument is ignored.}
}
\value{
  A results file chosen from a menu.
}
\note{
  AME: made changes so that options are compatible with those in 
  \code{load.allResFiles}. Previously, trouble occurred when overwriting.
}
\seealso{
  \code{\link{menuMenu}}, \code{\link{importCol2}}
}
\keyword{IO}
\keyword{file}% __ONLY ONE__ keyword per line


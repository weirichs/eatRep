\name{checkLEs}
\alias{checkLEs}
\title{Checks compatibility of linking errors with GADS data bases.}
\description{This function checks if a linking error \code{data.frame} is compatible with multiple trend \code{eatGADS} data bases. }
\usage{
checkLEs(filePaths, leDF)
}
\arguments{
  \item{filePaths}{
Character vectors with at least two paths to the \code{eatGADS} db files.
}
  \item{leDF}{
Linking error \code{data.frame}.
}
}
\details{
This function inspects whether all linking error variables correspond to variables in the \code{eatGADS} data base
and if the key variables also correspond to existing variables in the trend \code{eatGADS} data bases.
}
\value{
Returns a report list.
}

\examples{
trenddat1 <- system.file("extdata", "trend_gads_2010.db", package = "eatGADS")
trenddat2 <- system.file("extdata", "trend_gads_2015.db", package = "eatGADS")
trenddat3 <- system.file("extdata", "trend_gads_2020.db", package = "eatGADS")
load(system.file("extdata", "linking_error.rda", package = "eatRep"))
check1 <- checkLEs(c(trenddat1, trenddat2, trenddat3), lErr)
check2 <- checkLEs(c(trenddat1, trenddat2, trenddat3), lErr[1:14,])
}


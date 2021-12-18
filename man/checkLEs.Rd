\name{checkLEs}
\alias{checkLEs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Checks compatibility of linking errors with GADS data bases.}
\description{This function checks if a linking error \code{data.frame} is compatible with multiple trend \code{eatGADS} data bases. }
\usage{
checkLEs(filePaths, leDF)
%- maybe also 'usage' for other objects documented here.
}
\arguments{
  \item{filePaths}{
%%     ~~Describe \code{file} here~~
Character vectors with at least two paths to the \code{eatGADS} db files.
}
  \item{leDF}{
%%     ~~Describe \code{dif.term} here~~
Linking error \code{data.frame}. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function inspects whether all linking error variables correspond to variables in the \code{eatGADS} data base 
and if the key variables also correspond to existing variables in the trend \code{eatGADS} data bases.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a report list.
}

\examples{
# define eatGADs data bases
trenddat1 <- system.file("extdata", "trend_gads_2010.db", package = "eatGADS")
trenddat2 <- system.file("extdata", "trend_gads_2015.db", package = "eatGADS")
trenddat3 <- system.file("extdata", "trend_gads_2020.db", package = "eatGADS")
# use template linking Error Object
load(system.file("extdata", "linking_error.rda", package = "eatRep"))
check1 <- checkLEs(c(trenddat1, trenddat2, trenddat3), lErr)
check2 <- checkLEs(c(trenddat1, trenddat2, trenddat3), lErr[1:14,])
}


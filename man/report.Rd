\name{report}
\alias{report}

\title{Reporting function for all jk2-functions.}

\description{
Summarizes the output of the jk2-functions \code{jk2.mean}, \code{jk2.table},
\code{jk2.quantile}, and \code{jk2.glm}.
}

\usage{
report(jk2.out, trendDiffs = FALSE, add = list(),
      exclude = c("Ncases", "NcasesValid", "var"), printGlm = FALSE, digits = 3,
      printDeviance = FALSE)
}
\arguments{
  \item{jk2.out}{ output of the jk2-functions }
  \item{trendDiffs}{Logical: compute differences of trends? }
  \item{add}{Optional: additional columns for output. See examples of the jk2-functions}
  \item{exclude}{Which parameters should be excluded from reporting?}
  \item{printGlm}{Only relevant for \code{jk2.glm}: print summary on console?}
  \item{digits}{Only relevant if \code{printGlm} is TRUE. How many digits should be printed?}
  \item{printDeviance}{Only relevant if \code{printGlm} is TRUE. Should the deviance information printed additionally?}
}

\author{
Sebastian Weirich
}
\examples{
### see examples of the jk2-functions.
}

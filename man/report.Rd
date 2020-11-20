\name{report}
\alias{report}

\title{Reporting function for all jk2-functions.}

\description{
Summarizes the output of the jk2-functions \code{\link{repMean}},\code{\link{repTable}},
\code{\link{repQuantile}}, and \code{\link{repGlm}}.
}

\usage{
report(jk2.out, trendDiffs = FALSE, add = list(),
      exclude = c("Ncases", "NcasesValid", "var", "sampleSize"), printGlm = FALSE, round = TRUE,
      digits = 3, printDeviance = FALSE)
}
\arguments{
  \item{jk2.out}{ output of on eof the four rep-functions. }
  \item{trendDiffs}{Logical: compute differences of trends? }
  \item{add}{Optional: additional columns for output. See examples of the jk2-functions}
  \item{exclude}{Which parameters should be excluded from reporting?}
  \item{printGlm}{Only relevant for \code{\link{repGlm}}: print summary on console?}
  \item{round}{Logical: should the results be rounded to a limited number of digits?}
  \item{digits}{How many digits should be used for rounding?}
  \item{printDeviance}{Only relevant for \code{\link{repGlm}} and if \code{printGlm} is TRUE.
  Should the deviance information printed additionally? Note: To print deviance information,
  the argument \code{poolMethod} of the \code{\link{repGlm}} function must be set to \code{"scalar"}. }
}

\author{
Benjamin Becker, Sebastian Weirich
}
\examples{
### see examples of the jk2-functions.
}

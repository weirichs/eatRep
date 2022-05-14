\name{report}
\alias{report}

\title{Reporting function for \code{\link{repMean}},\code{\link{repTable}},
\code{\link{repQuantile}}, and \code{\link{repGlm}} }

\description{
Summarizes the output of the four main functions \code{\link{repMean}},\code{\link{repTable}},
\code{\link{repQuantile}}, and \code{\link{repGlm}}, and provides a single data.frame with all
results.
}

\usage{
report(repFunOut, trendDiffs = FALSE, add = list(),
      exclude = c("Ncases", "NcasesValid", "var", "sampleSize"), printGlm = FALSE,
      round = TRUE, digits = 3, printDeviance = FALSE,
      target = c("default", "BT2021", "BT2021.table"),
      wholeGroupName = "Deutschland")
}
\arguments{
  \item{repFunOut}{output of one of the four \code{eatRep}-functions. }
  \item{trendDiffs}{Logical: compute differences of trends? }
  \item{add}{Optional: additional columns for output. See examples of the jk2-functions}
  \item{exclude}{Which parameters should be excluded from reporting?}
  \item{printGlm}{Only relevant for \code{\link{repGlm}}: print summary on console?}
  \item{round}{Logical: should the results be rounded to a limited number of digits?}
  \item{digits}{How many digits should be used for rounding?}
  \item{printDeviance}{Only relevant for \code{\link{repGlm}} when other than the
  identical function is used as link function, and if \code{printGlm} is TRUE.
  Should the deviance information printed additionally? Note: To print deviance information,
  the argument \code{poolMethod} of the \code{\link{repGlm}} function must be set
  to \code{"scalar"}. }
  \item{target}{The target format for the reported data.frame. Option \code{"BT2021"}
  optimizes the \code{\link{repMean}} and \code{\link{repTable}} output format according
  to the BT2021 study. Option \code{"BT2021.table"} additionally adapts the
  \code{\link{repTable}} output format according to the BT2021 study. For general purposes,
  we recommend to use the \code{"default"} option.}
  \item{wholeGroupName}{Optional: only relevant if \code{target} is \code{"BT2021"} or \code{"BT2021.table"}.
  Which notation should be used for whole population in the output?}
}
\value{
A data frame with at least nine columns.
\item{group}{Denotes the group an analysis belongs to. If no groups were specified and/or analysis for the
whole sample were requested, the value of \sQuote{group} is \sQuote{wholeGroup}.}
\item{depVar}{Denotes the name of the dependent variable in the analysis. }
\item{modus}{Denotes the mode of the analysis. For example, if a JK2 regression analysis was
conducted, \sQuote{modus} takes the value \sQuote{JK2.glm}. If a mean analysis without any replicates
was conducted, \sQuote{modus} takes the value \sQuote{CONV.mean}.}
\item{comparison}{Denotes whether group mean comparisons or cross-level comparisons were conducted.
Without any comparisons, \sQuote{comparison} takes the value \sQuote{NA} }
\item{parameter}{Denotes the parameter of the corresponding analysis. If regression analysis was applied,
the regression parameter is given. Amongst others, the \sQuote{parameter} column takes the values
\sQuote{(Intercept)} and \sQuote{gendermale} if \sQuote{gender} was the independent variable, for instance.
If mean analysis was applied, the \sQuote{parameter} column takes the values \sQuote{mean}, \sQuote{sd},
\sQuote{var}, or \sQuote{Nvalid}. See the examples of \code{\link{repMean}},\code{\link{repTable}},
\code{\link{repQuantile}}, or \code{\link{repGlm}} for further details. }
\item{depVar}{Denotes the name of the dependent variable (only if \code{\link{repGlm}} was called before) }
\item{est}{Denotes the estimate of the corresponding analysis.}
\item{se}{Denotes the standard error of the corresponding estimate.}
\item{p}{Denotes the p value of the estimate.}
}
\author{
Benjamin Becker, Sebastian Weirich
}
\examples{
### see examples of the eatRep main functions.
}

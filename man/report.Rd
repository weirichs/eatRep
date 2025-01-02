\name{report}
\alias{report}
\alias{report2}

\title{Reporting functions for \code{\link{repMean}}, \code{\link{repTable}},
\code{\link{repQuantile}}, and \code{\link{repGlm}} }

\description{
Summarizes the output of the four main functions \code{\link{repMean}}, \code{\link{repTable}},
\code{\link{repQuantile}}, and \code{\link{repGlm}}, and provides a single data.frame with all
results.}

\usage{
report(repFunOut, trendDiffs = deprecated(), add=list(),
       exclude = c("NcasesValid", "var"),
       printGlm = FALSE, round = TRUE, digits = 3, printDeviance = FALSE,
       printSE_correction = FALSE)
report2(repFunOut, add=list(), exclude = c("NcasesValid", "var"), printGlm = FALSE,
       round = TRUE, digits = 3, printDeviance = FALSE, printSE_correction = FALSE) }
\arguments{
  \item{repFunOut}{output of one of the four \code{eatRep}-functions. }
  \item{trendDiffs}{deprecated. In earlier versions of the package, this argument was used
  to determine differences in trends. As differences in trends are equivalent to the trend
  of differences (no matter whether group or cross-level differences), the argument was
  deprecated. If the user specifies group or cross-level difference along with trends,
  trends of differences are computed as well.}
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
  \item{printSE_correction}{Logical: Print the differences of original SEs of cross
  differences (method \code{"old"}) and SEs obtained by the \code{"wec"} or \code{"rep"} method.}
}
\value{
\code{report} and \code{report2} differ in the output which is returned. The output of \code{report2}
is optimized for further processing, i.e. drawing plots by means of the \code{eatPlot} package. For
\code{report}, the output is a data frame with at least nine columns.
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
For \code{report2}, the output is a list with four data.frames. The first data.frame \code{plain} summarizes the
results in human-readable form. The data.frames 2 to 4 (\code{comparisons}, \code{group}, \code{estimate}) are redundant
to \code{plain} and contain the results in a technical presentation suitable for further processing in \code{eatPlot}.
\item{plain}{The complete results in human-readable form.}
\item{comparison}{An allocation table that indicates which comparison (group comparison or cross-level comparison)
relates to which groups.}
\item{group}{A table that assigns an ID to each analysis unit. This makes it easier later on to read from the output
which comparison relates to which groups. This simplifies the assignment, especially when comparing comparisons
(i.e., cross-level differences of group differences).}
\item{estimate}{The results of the analyses, assigned to their IDs.}
}
\author{
Benjamin Becker, Sebastian Weirich
}
\examples{
### see examples of the eatRep main functions.
}

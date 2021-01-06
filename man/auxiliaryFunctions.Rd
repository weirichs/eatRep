\name{auxiliaryFunctions}
\alias{funAdjust}
\alias{funAdjustEL}
\alias{groupVersusTotalMean}
\alias{getOutputIfSingular}
\alias{getOutputIfSingularT1}
\alias{funadjustLavaanWec}
\alias{getOutputIfSingularWec}
\title{
Internally used auxiliary functions
}
\description{
These functions are used internally and have to be exported to the namespace
to be callable by \code{\link[survey]{withReplicates}}, for example. The user
should not call these functions directly.
}
\usage{
funAdjust (d, x, a, w)
funAdjustEL(d, x, a, w)
groupVersusTotalMean(d, g, w)
getOutputIfSingularWec ( glmRes )
getOutputIfSingular ( glmRes )
getOutputIfSingularT1 ( glmRes)
funadjustLavaanWec (d, x, i, w)
}
\arguments{
\item{d}{dependent variable, provided as vector}
\item{x}{grouping variable, provided as vector}
\item{a}{matrix of variables used for adjusted means}
\item{w}{optional: weighting variable, provided as vector}
\item{g}{grouping variable, provided as vector}
\item{i}{independent variable, provided as vector}
\item{glmRes}{results object of the \code{glm} function}
}
\value{
data frame with two columns: estimates and standard errors, both only used for further
processing in \code{eatRep}.
}

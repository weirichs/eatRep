\name{pool.R2}
\alias{pool.R2}
\title{Compute \eqn{R^2} in multiple imputed and nested multiple imputed data}
\description{With (nested) multiple imputations, the determination coefficient \eqn{R^2} has to be
computed for each imputed data set and pooled afterwards. \code{pool.R2} provide pooling routines
according to Harel (2009). The function requires that the \eqn{R^2} coefficients from the multiple imputed
analyses are already available.}
\usage{
pool.R2 ( r2, N, verbose = TRUE )
}
\arguments{
  \item{r2}{
For multiple imputed data, a numeric vector of \eqn{R^2} values. For nested multiple imputed data, a
list of numeric vectors of \eqn{R^2} values. The number of list elements must correspond to the number of nests.
The number of \eqn{R^2} values within each list element must be equal and must correspond to the number of
imputations within each nest.
}
  \item{N}{
Optional: the sample size of each imputed data set. Only necessary if the standard error for the pooled \eqn{R^2}
should be computed. The structure of the \code{N} object must correspond to the structure of the \code{r2}
object. See examples for further details.
}
  \item{verbose}{
Optional: Print additional messages to console?
}
}
\value{
Returns a data.frame with one or two columns which contains the pooled \eqn{R^2} value and optionally it's standard error.
}
\references{
  Harel, O. (2009): The estimation of \eqn{R^2} and adjusted \eqn{R^2} in incomplete data
  sets using multiple imputation. \emph{Journal of Applied Statistics.} \bold{36, 10}, 1109--1118.
}
\examples{
# multiple imputation, assume that the regression analysis was fitted for five imputed data sets,
# resulting in five R^2 values. Assume sample sizes of 340
r2 <- c(0.12395, 0.15261, 0.16125, 0.11029, 0.1871)
Ns <- rep(340,5)
pool.R2 ( r2=r2, N=Ns)
# without standard error
pool.R2 ( r2=r2)
# nested multiple imputation
r2 <- list(nest1 = c(0.12395, 0.15261, 0.16125, 0.11029, 0.1871),
           nest2 = c(0.10603, 0.08876, 0.09248, 0.13331, 0.1114),
           nest3 = c(0.17228, 0.25203, 0.13132, 0.23331, 0.10069))
Ns <- lapply(1:3, FUN = function (x) {rep(290, 5)})
pool.R2 ( r2=r2, N=Ns)
# without standard error
pool.R2 ( r2=r2)
}

\name{eatRep-package}
\alias{eatRep-package}
\docType{package}
\title{
	Statistical analyses in complex survey designs with multiple imputed data and trend estimation.
}
\description{
  The package provide functions to computes some basic statistic operations---(adjusted) means, standard deviations,
  frequency tables, percentiles and generalized linear models---in complex survey designs comprising multiple
  imputed variables and/or a clustered sampling structure which both deserve special procedures at least in
  estimating standard errors. In large-scale assessments, standard errors are comprised of three components:
  the measurement error, the sampling error, and (if trend estimation of at least two times of measurement
  are involved) the linking error.

  \strong{Measurement error:} In complex surveys or large-scale assessments, measurement errors are taken
  into account by the mean of multiple imputed variables. The computation of standard errors for the mean
  of a multiple imputed variable (e.g. plausible values) involves the formulas provided by Rubin (1987).
  Computing standard errors for the mean of a nested imputed variable involves the formulas provided by
  Rubin (2003). Both methods are implemented in the package. The estimation of \eqn{R^2} and adjusted
  \eqn{R^2} in linear and generalized linear regression models with multiple imputed data sets is
  realized using the methods provided in Harel (2009).

  \strong{Sampling error:} Computation of sampling errors of variables which stem from a clustered design may
  involve replication methods like balanced repeated replicate (BRR), bootstrap or Jackknife methods.
  See Westat (2000), Foy, Galia & Li (2008), Rust and Rao (1996), and Wolter (1985) for details. To date,
  the Jackknife-1 (JK1), Jackknife-2 (JK2) and the Balanced Repeated Replicates (BRR; optionally with Fay's
  method) procedures are supported.

  \strong{Linking error:} Lastly, standard errors for trend estimates may involve incorporating
  linking errors to account for potential differential item functioning or item parameter drift.
  \code{eatRep} allows to account for linking error when computing standard errors for trend
  estimates. Standard error estimation is conducted according to the operational practice in
  PISA, see equation 5 in Sachse & Haag (2017).

  The package \code{eatRep} is designed to combine one or several error types which is necessary,
  for example, if (nested) multiple imputed data are used in clustered designs. Considering the
  structure is relevant especially for the estimation of standard errors. The estimation of national
  trends requires a sequential analysis for both measurements and a comparison of estimates between them.

  Technically, \code{eatRep} is a wrapper for the \code{survey} package (Lumley, 2004). Each function in
  \code{eatRep} corresponds to a specific function in \code{survey} which is called repeatedly during the analysis.
  Hence, a nested loop is used. We use \dQuote{trend replicates} in the outer loop, \dQuote{imputation replicates} 
  in the middle loop to account for multiple imputed data, and \dQuote{cluster replicates} in the inner loop to 
  account for the clustered sampling structure. While the functional principle of \code{survey} is based on 
  replication of standard analyses, \code{eatRep} is based on replication of \code{survey} analyses to take 
  multiple imputed data into account. More recent versions of the package additionally allow estimations using
  the \code{BIFIEsurvey} package instead of \code{survey} which provide substantial advantages in terms of speed. 
  
  For each imputed data set in each measurement, i.e. in the inner loop, the \code{eatRep} function first creates 
  replicate weights based on the primary sampling unit (PSU) variable and the replication indicator variable. In 
  the jackknife procedure, the first one is often referred to as \dQuote{jackknife zone}, whereas the second one 
  is often referred to as \dQuote{jackknife replicate}. The number of distinct units in the PSU variable defines
  the number of replications which are necessary due to the clustered structure. A design object is created and 
  the appropriate \code{survey} function is called. The process is repeated for each imputed dataset and the 
  results of the analyses are pooled. The pooling procedure varies in relation to the type of variable to be 
  pooled. For examples, means or regression coefficients are pooled according to Rubin (1987) or Rubin (2003). 
  \eqn{R^2} is pooled according to Harel (2009), using a Fisher \emph{z}-transformation. Chi-square distributed values 
  are pooled according to Thomas and Rao (1990) for clustered data and according to Enders (2010) and 
  Allison (2002) for multiple imputed data. For trend analyses, the whole process is repeated two times 
  (according to the two measurements) and the difference of the estimates are computed along with their 
  pooled standard errors. 
  
  Without trend estimation, the outer loop has only one cycle (instead of two). Without multiple imputations, 
  the middle loop has only one cycle. Without a clustered sampling structure (i.e, in a random sample), the 
  inner loop has only one cycle. Without trend, imputation and clustered structure, no replication is performed 
  at all. To compute simple mean estimates, for example, \code{eatRep} then simply calls \code{mean} instead 
  of \code{svymean} from the \code{survey} package. A special case occurs with nested multiple imputation. 
  We then have four loops in a nested structure. Hence, the corresponding analyses may take considerably 
  computational effort. 
 
  \emph{Important note:} Starting with version 0.10.0, several methods for the standard error estimation
  of cross level differences are implemented. Prior to version 0.10.0, the standard error for the difference
  between one single group (e.g., Belgium) and the total population (which is comprised of several states including 
  Belgium) was estimated as if both groups would have been independent from each other. The standard errors,
  however, are biased then. Two new methods are now applicable using the argument \code{crossDiffSE} in 
  \code{\link{repMean}} and provide unbiased standard errors---weighted effect coding (wec) and replication
  methods (rep); see, for example te Grotenhuis et al. (2017) and Weirich et al. (2021). The old method is still available by
  using \code{crossDiffSE = "old"}. Note that the default method now is weighted effect coding.

  \emph{Second important note:} Starting with version 0.13.0, function names have been changed due to
  inconsistent former denomination: Function \code{jk2.mean} now goes under the name of \code{\link{repMean}},
  \code{jk2.table} was  renamed to \code{\link{repTable}}, \code{jk2.quantile} was  renamed to \code{\link{repQuantile}},
  and \code{jk2.glm} now goes under the name of \code{\link{repGlm}}. The old functions are deprecated and will
  be removed in further package publications. Renaming was driven by the fact that the corresponding
  functions now have broader range of methods than only jackknife-2.
}
\details{
\tabular{ll}{
Package: \tab eatRep\cr
Type: \tab Package\cr
Version: \tab 0.14.5\cr
Date: \tab 2022-07-05\cr
License: \tab GPL(>=2)
}
}
\author{
    Authors: Sebastian Weirich <sebastian.weirich@iqb.hu-berlin.de>, Martin Hecht <martin.hecht@hu-berlin.de>,
    Benjamin Becker <b.becker@iqb.hu-berlin.de>
}
\references{
  Allison, P. D. (2002). Missing data. Newbury Park, CA: Sage.

  Enders, C. K. (2010). Applied missing data analysis. Guilford Press.

  Foy, P., Galia , J. & Li, I. (2008). Scaling the data from the TIMSS 2007 mathematics
  and science assessment. In J. F. Olson, M. O. Martin & I. V. S. Mullis (ed.),
  \emph{TIMSS 2007 Technical Report} (S. 225--280). Chestnut Hill, MA: TIMSS & PIRLS
  International Study Center, Lynch School of Education, Boston College.
  
  Harel, O. (2009): The estimation of \eqn{R^2} and adjusted \eqn{R^2} in incomplete data 
  sets using multiple imputation. \emph{Journal of Applied Statistics.} \bold{36, 10}, 1109--1118.

  Lumley, T. (2004). Analysis of complex survey samples. \emph{Journal of Statistical Software} \bold{9(1)}: 1--19

  Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys.} New York: Wiley.

  Rubin, D.B. (2003): Nested multiple imputation of NMES via partially incompatible MCMC.
  \emph{Statistica Neerlandica} \bold{57, 1}, 3--18.
  
  Rust, K., & Rao, JNK. (1996): Variance estimation for complex surveys using 
  replication techniques. \emph{Statistical Methods in Medical Research} \bold{5}, 283--310.

  Sachse, K. A. & Haag, N. (2017). Standard errors for national trends in international
  large-scale assessments in the case of cross-national differential item functioning. \emph{Applied
  Measurement in Education, 30}, (2), 102-116. http://dx.doi.org/10.1080/08957347.2017.1283315

	Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics
		and standard errors in covariance structure analysis.

  te Grotenhuis, M., Pelzer, B., Eisinga, R., Nieuwenhuis, R., Schmidt-Catran, A., & Konig, R. (2017).
  When size matters: advantages of weighted effect coding in observational studies.
  \emph{International Journal of Public Health.} \bold{62}, 163--167.

  Thomas, D. R. & Rao, JNK (1990): Small-sample comparison of level and power for simple goodness-of-
  fit statistics under cluster sampling. JASA 82:630-636
  
  Weirich, S., Hecht, M., Becker, B. et al. (2021). Comparing group means with the total mean in random samples,
  surveys, and large-scale assessments: A tutorial and software illustration. Behavior Research Methods.
  https://doi.org/10.3758/s13428-021-01553-1

  Westat (2000). \emph{WesVar.} Rockville, MD: Westat.

  Wolter, K. M. (1985). \emph{Introduction to variance estimation.} New York: Springer.
}
\keyword{ package }

\name{eatRep-package}
\alias{eatRep-package}
\docType{package}
\title{
	Statistical analyses in complex survey designs with multiple imputed data and trend estimation.
}
\description{
  Computes some basic statistic operations (means, standard deviations, frequency tables,
  percentiles and generalized linear models) in complex survey designs comprising multiple imputed variables
  and/or a clustered sampling structure which both deserve special procedures at least in estimating standard errors.

  For example, computing standard errors for the mean of a multiple imputed variable (e.g. plausible values) involves the
  formulas provided by Rubin (1987). Computing standard errors for the mean of a nested imputed variable
  involves the formulas provided by Rubin (2003). Both methods are implemented in the package. The estimation of 
  \eqn{R^2} and adjusted \eqn{R^2} in linear and generalized linear regression models with multiple imputed data sets is 
  realized using the methods provided in Harel (2009). 

  Moreover, computing standard errors for the mean of a variable which stems from a clustered design may involve
  replication methods like balanced repeated replicate (BRR), bootstrap or Jackknife methods.
  See Weststat (2000), Foy, Galia & Li (2008), Rust and Rao (1996), and Wolter (1985) for details. 
  To date, the Jackknife-1 (JK1), Jackknife-2 (JK2) and the Balanced Repeated Replicates (BRR) procedures are supported.

  The package \code{eatRep} is designed to combine both methods which is necessary if (nested) multiple imputed
  data are used in clustered designs. Considering the structure is relevant especially for the estimation of
  standard errors. The estimation of national trends requires a sequential analysis for both measurements
  and a comparison of estimates between them. 

  Technically, \code{eatRep} is a wrapper for the \code{survey} package (Lumley, 2004). Each function in
  \code{eatRep} corresponds to a specific function in \code{survey} which is called repeatedly during the analysis.
  Hence, a nested loop is used. We use \dQuote{trend replicates} in the outer loop, \dQuote{imputation replicates} 
  in the middle loop to account for multiple imputed data, and \dQuote{cluster replicates} in the inner loop to 
  account for the clustered sampling structure. While the functional principle of \code{survey} is based on 
  replication of standard analyses, \code{eatRep} is based on replication of \code{survey} analyses to take 
  multiple imputed data into account. 
  
  For each imputed data set in each measurement, i.e. in the inner loop, the \code{eatRep} function first creates 
  replicate weights based on the primary sampling unit (PSU) variable and the replication indicator variable. In 
  the jackknife procedure, the first one is often referred to as \dQuote{Jackknife Zone}, whereas the second one 
  is often referred to as \dQuote{Jackknife Replicate}. The number of distinct units in the PSU variable define 
  the number of replications which are necessary due to the clustered structure. A design object is created and 
  the appropriate \code{survey} function is called. The process is repeated for each imputed dataset and the 
  results of the analyses are pooled. The pooling procedure varies in relation to the type of variable to be 
  pooled. For examples, means or regression coefficients are pooled according to Rubin (1987) or Rubin (2003). 
  \eqn{R^2} is pooled according to Harel (2009), using a Fisher z-transformation. Chi-square distributed values 
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
 
  \emph{Important note:} The structure of the the \code{eatRep} functions varied substantially between versions
  0.5.0 and 0.6.0. Up to version 0.5.0, the data has to be provided in the wide format. Beginning with version 
  0.6.0, \code{eatRep} functions need the long format. This distinction practically means that version 0.5.0
  allows to analyze data where, for example, the number of imputations is different between independent and 
  dependent variables, albeit the second one is \emph{not} nested within the first one. This case is conceptually 
  questionable and it is not clear how to imply the pooling rules. Hence, this is no longer supported in version 
  0.6.0 and higher. The number of imputations have to be equal or a nested structure must be guaranteed. 
}
\details{
\tabular{ll}{
Package: \tab eatRep\cr
Type: \tab Package\cr
Version: \tab 0.9.4\cr
Date: \tab 2018-10-04\cr
License: \tab GPL(>=2)
}
}
\author{
    Author/maintainer: Sebastian Weirich <sebastian.weirich@iqb.hu-berlin.de>
}
\references{
  Allison, P. D. (2002). Missing data. Newbury Park, CA: Sage.

  Enders, C. K. (2010). Applied missing data analysis. Guilford Press.

  Foy, P., Galia , J. & Li, I. (2008). Scaling the data from the TIMSS 2007 mathematics
  and science asssessment. In J. F. Olson, M. O. Martin & I. V. S. Mullis (ed.),
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

	Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics
		and standard errors in covariance structure analysis.

  Thomas, D. R. & Rao, JNK (1990): Small-sample comparison of level and power for simple goodnessof-
  fit statistics under cluster sampling. JASA 82:630-636

  Westat (2000). \emph{WesVar.} Rockville, MD: Westat.

  Wolter, K. M. (1985). \emph{Introduction to variance estimation.} New York: Springer.

}
\keyword{ package }
\seealso{
}
\examples{
}

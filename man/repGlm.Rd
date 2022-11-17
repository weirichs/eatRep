\name{repGlm}
\alias{repGlm}
\alias{jk2.glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Replication methods (JK1, JK2 and BRR) for linear regression models and trend estimation.}
\description{Compute generalized linear models for complex cluster designs with multiple imputed variables based 
on the Jackknife (JK1, JK2) or balanced repeated replicates (BRR) procedure. Conceptually, the function combines replication 
methods and methods for multiple imputed data. Technically, this is a wrapper for the \code{\link[survey]{svyglm}} function
of the \code{survey} package.
}
\usage{
repGlm(datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"), PSU = NULL,
        repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL,
        group.splits = length(groups), group.delimiter = "_",
        cross.differences = FALSE, trend = NULL, linkErr = NULL, formula,
        family=gaussian, forceSingularityTreatment = FALSE,
        glmTransformation = c("none", "sdY"), doCheck = TRUE, na.rm = FALSE,
        poolMethod = c("mice", "scalar"), useWec = FALSE,
        scale = 1, rscales = 1, mse=TRUE, rho=NULL, hetero=TRUE,
        se_type = c("HC3", "HC0", "HC1", "HC2", "CR0", "CR2"),
        clusters = NULL, crossDiffSE.engine= c("lavaan", "lm"),
        stochasticGroupSizes = FALSE, verbose = TRUE, progress = TRUE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datL}{
%%     ~~Describe \code{file} here~~
Data frame in the long format (i.e. each line represents one ID unit in one imputation of one nest) containing all 
variables for analysis.
}
  \item{ID}{
%%     ~~Describe \code{dif.term} here~~
Variable name or column number of student identifier (ID) variable. ID variable must not contain any missing values. 
}
  \item{wgt}{
%%     ~~Describe \code{split.dif} here~~
Optional: Variable name or column number of weighting variable. If no weighting variable is specified,
all cases will be equally weighted.
}
  \item{type}{
%%     ~~Describe \code{abs.dif.bound} here~~
Defines the replication method for cluster replicates which is to be applied. Depending on \code{type}, additional
arguments must be specified (e.g., \code{PSU} and/or \code{repInd} or \code{repWgt}).
}
  \item{PSU}{
%%     ~~Describe \code{abs.dif.bound} here~~
Variable name or column number of variable indicating the primary sampling unit (PSU). When a jackknife procedure is applied, 
the PSU is the jackknife zone variable. If \code{NULL}, no cluster structure is assumed and
standard errors are computed according to a random sample.
}
  \item{repInd}{
%%     ~~Describe \code{sig.dif.bound} here~~
Variable name or column number of variable indicating replicate ID. In a jackknife procedure, this is the jackknife replicate 
variable. If \code{NULL}, no cluster structure is assumed and standard errors are computed according to a random sample.
}
  \item{repWgt}{
%%     ~~Describe \code{sig.dif.bound} here~~
Normally, replicate weights are created by \code{repGlm} directly from \code{PSU} and \code{repInd} variables. Alternatively,
if replicate weights are included in the data.frame, specify the variable names or column number in the \code{repWgt} argument.
}
\item{nest}{Optional: name or column number of the nesting variable. Only applies in nested multiple imputed data sets.
}
\item{imp}{Optional: name or column number of the imputation variable. Only applies in multiple imputed data sets.
}
  \item{groups}{
%%     ~~Describe \code{sig.dif.bound} here~~
Optional: vector of names or column numbers of one or more grouping variables. 
}
  \item{group.splits}{
%%     ~~Describe \code{sig.dif.bound} here~~
Optional: If groups are defined, \code{group.splits} optionally specifies whether analysis should be done also
in the whole group or overlying groups. See examples for more details.
}
  \item{group.delimiter}{
%%     ~~Describe \code{sig.dif.bound} here~~
Character string which separates the group names in the output frame.
}
  \item{cross.differences}{
%%     ~~Describe \code{sig.dif.bound} here~~
Either a list of vectors, specifying the pairs of levels for which cross-level differences should be computed.
Alternatively, if \code{TRUE}, cross-level differences for all pairs of levels are computed. If \code{FALSE}, no cross-level
differences are computed. (see examples 2a, 3, and 4 in the help file of the \code{repMean} function)
}
  \item{trend}{
%%     ~~Describe \code{sig.dif.bound} here~~
Optional: name or column number of the trend variable. Note: Trend variable must have exact two levels. Levels for 
grouping variables must be equal in both 'sub populations' partitioned by the trend variable. 
}
  \item{linkErr}{
%%     ~~Describe \code{sig.dif.bound} here~~
Optional: name or column number of the linking error variable. If \code{NULL}, a linking error of 0 will be assumed in trend estimation.
}
  \item{formula}{
%%     ~~Describe \code{sig.dif.bound} here~~
Model formula, see help page of \code{glm} for details. 
}
  \item{family}{
%%     ~~Describe \code{sig.dif.bound} here~~
A description of the error distribution and link function to be used in the model. See help page of \code{glm} for details. 
}
  \item{forceSingularityTreatment}{
%%     ~~Describe \code{sig.dif.bound} here~~
Logical: Forces the function to use the workaround to handle singularities in regression models.
}
  \item{glmTransformation}{
Optional: Allows for transformation of parameters from linear regression and logistic regression before pooling.
Useful to compare parameters from different glm models, see Mood (2010). Note: This argument applies only if 
\code{forceSingularityTreatment} is set to 'TRUE'.
}
  \item{doCheck}{
Logical: Check the data for consistency before analysis? If \code{TRUE} groups with insufficient data are excluded from
analysis to prevent subsequent functions from crashing.
}
  \item{na.rm}{
Logical: Should cases with missing values be dropped?
}
  \item{poolMethod}{
Which pooling method should be used? The ``mice'' method is recommended.
}
  \item{useWec}{
Logical: use weighted effect coding? 
}
  \item{scale}{
scaling constant for variance, for details, see help page of \code{svrepdesign} from the \code{survey} package
}
  \item{rscales}{
scaling constant for variance, for details, see help page of \code{svrepdesign} from the \code{survey} package
}
  \item{mse}{
Logical: If \code{TRUE}, compute variances based on sum of squares around the point estimate, rather than the mean of the replicates.
See help page of \code{\link[survey]{svrepdesign}} from the \code{survey} package for further details.
}
  \item{rho}{
Shrinkage factor for weights in Fay's method. See help page of \code{\link[survey]{svrepdesign}} from the \code{survey} package for further details.
}
  \item{hetero}{
Logical: Assume heteroscedastic variance for weighted effect coding? Only applies for random samples, i.e. if no replication analyses are executed. 
}
  \item{se_type}{
The sort of standard error sought for cross level differences. Only applies if \code{crossDiffSE == "wec"} and \code{hetero == TRUE}
and \code{crossDiffSE.engine == "lm"}. See the help page of \code{\link[estimatr]{lm_robust}} from the \code{estimatr} package for further details.
}
  \item{clusters}{
Optional: Variable name or column number of cluster variable. Only necessary if weighted effecting coding
should be performed using heteroscedastic variances. See the help page of \code{\link[estimatr]{lm_robust}}
from the \code{estimatr} package for further details.
}
  \item{crossDiffSE.engine}{
Optional: Sort of estimator which should be used for standard error estimation in weighted effect coding regression.
Only applies if \code{useWec == TRUE}. To date, only lavaan allows for stochastic group sizes.
}
  \item{stochasticGroupSizes}{
Logical: Assume stochastic group sizes for using weighted effect coding regression with categorical predictors? Note: To date,
only lavaan allows for stochastic group sizes. Stochastic group sizes cannot be assumed if any replication method
(jackknife, BRR) is applied.
}
  \item{verbose}{
Logical: Show analysis information on console?
}
  \item{progress}{
Logical: Show progress bar on console?
}
}
\details{
Function first creates replicate weights based on \code{PSU} and \code{repInd} variables according to JK2 or
BRR procedure. According to multiple imputed data sets, a workbook with several analyses is created.
The function afterwards serves as a wrapper for \code{\link[survey]{svyglm}} implemented in the \code{survey} package.
The results of the several analyses are then pooled according to Rubin's rule, which is adapted for nested 
imputations if the \code{nest} argument implies a nested structure.
}
\value{
A list of data frames in the long format. The output can be summarized using the \code{report} function.
The first element of the list is a list with either one (no trend analyses) or two (trend analyses)
data frames with at least six columns each. For each subpopulation denoted by the \code{groups}
statement, each dependent variable, each parameter and each coefficient the corresponding value is given. 
\item{group}{Denotes the group an analysis belongs to. If no groups were specified and/or analysis for the 
whole sample were requested, the value of \sQuote{group} is \sQuote{wholeGroup}.}
\item{depVar}{Denotes the name of the dependent variable in the analysis. }
\item{modus}{Denotes the mode of the analysis. For example, if a JK2 analysis without sampling weights was 
conducted, \sQuote{modus} takes the value \sQuote{jk2.unweighted}. If a analysis without any replicates but with sampling
weights was conducted, \sQuote{modus} takes the value \sQuote{weighted}.}
\item{parameter}{Denotes the parameter of the regression model for which the corresponding value is given
further. Amongst others, the \sQuote{parameter} column takes the values \sQuote{(Intercept)} and \sQuote{gendermale} if \sQuote{gender} 
was the dependent variable, for instance. See example 1 for further details.}
\item{coefficient}{Denotes the coefficient for which the corresponding value is given further. Takes the 
values \sQuote{est} (estimate) and \sQuote{se} (standard error of the estimate).}
\item{value}{The value of the parameter estimate in the corresponding group.}
If groups were specified, further columns which are denoted by the group names are added to the data frame. 
}
\references{
te Grotenhuis, M., Pelzer, B., Eisinga, R., Nieuwenhuis, R., Schmidt-Catran, A., & Konig, R. (2017).
When size matters: advantages of weighted effect coding in observational studies. \emph{International Journal of Public Health.} \bold{62}, 163--167.
}
\examples{
### load example data (long format)
data(lsa)
### use only the first nest
bt         <- lsa[which(lsa[,"nest"] == 1),]
### use only data from 2010
bt2010     <- bt[which(bt[,"year"] == 2010),]
## use only reading data
bt2010read <- bt2010[which(bt2010[,"domain"] == "reading"),]

### Example 1: Computes linear regression from reading score on gender separately
### for each country. Assume no nested structure.
mod1 <- repGlm(datL = bt2010read, ID = "idstud", wgt = "wgt", type = "jk2",
        PSU = "jkzone", repInd = "jkrep", imp = "imp", groups = "country",
        formula = score~sex, family ="gaussian")
res1 <- report(mod1, printGlm = TRUE)

\donttest{
### Example 2: Computes log linear regression from pass/fail on ses and gender
### separately for each country in a nested structure. Assuming equally weighted
### cases by omitting "wgt" argument
dat  <- lsa[intersect(which(lsa[,"year"] == 2010), which(lsa[,"domain"] == "reading")),]
mod2 <- repGlm(datL = dat, ID = "idstud", type = "JK2",  PSU = "jkzone",
        repInd = "jkrep", imp = "imp", nest="nest", groups = "country",
        formula = passReg~sex*ses, family = quasibinomial(link="logit"))
res2 <- report(mod2, printGlm = TRUE)

### Example 3: Like example 1, but without any replication methods
### trend estimation (without linking error) and nested imputation
dat  <- lsa[which(lsa[,"domain"] == "reading"),]
mod3 <- repGlm(datL = dat, ID = "idstud", wgt = "wgt", imp = "imp", nest = "nest",
        groups = "country",  formula = score~sex, trend = "year")
res3 <- report(mod3, printGlm = TRUE)

### Example 4: weighted effect coding to estimate whether a specific country's mean
### differs from the overall mean (whereas the overall population is a composite of 
### all countries). The procedure adapts the weighted effect coding procedures 
### described in te Grotenhuis (2017) for multiple imputation and replication methods. 
mod4 <- repGlm(datL = bt2010read, ID = "idstud", wgt = "wgt", type = "jk2",
        PSU = "jkzone", repInd = "jkrep", imp = "imp", formula = score~country,
        useWec=TRUE)
res4 <- report(mod4, printGlm = FALSE)
}}


\name{repMean}
\alias{repMean}
\alias{jk2.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Replication methods (JK1, JK2 and BRR) for descriptive statistics.}
\description{Compute totals, means, adjusted means, mean differences, variances and standard deviations
with standard errors in random or clustered or complex samples. Variance estimation in complex cluster
designs based on Jackknife (JK1, JK2) or Balanced Repeated Replicates (BRR) procedure. Moreover, analyses
can be customized for multiple or nested imputed variables, applying the combination rules of Rubin (1987)
for imputed data and Rubin (2003) for nested imputed data. Conceptually, the function combines replication
methods and methods for multiple imputed data. Trend estimation as usual in large-scale assessments is supported as well.
Technically, this is a wrapper for the \code{\link[survey]{svymean}} and \code{\link[survey]{svyvar}} functions of the \code{survey} package.}
\usage{
repMean (datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"), PSU = NULL,
         repInd = NULL, jkfac=NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL,
         group.splits = length(groups), group.differences.by = NULL,
         cross.differences = FALSE, crossDiffSE = c("wec", "rep","old"),
         adjust = NULL, useEffectLiteR = FALSE, nBoot = 100, group.delimiter = "_",
         trend = NULL, linkErr = NULL, dependent, na.rm = FALSE, doCheck = TRUE,
         engine = c("survey", "BIFIEsurvey"), scale = 1, rscales = 1, mse=TRUE,
         rho=NULL, hetero=TRUE, se_type = c("HC3", "HC0", "HC1", "HC2", "CR0", "CR2"),
         clusters = NULL, crossDiffSE.engine= c("lavaan", "lm"),
         stochasticGroupSizes = FALSE, verbose = TRUE, progress = TRUE, nCores=NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datL}{
Data frame in the long format (i.e. each line represents one ID unit in one imputation of one nest) containing all
variables for analysis.
}
  \item{ID}{
Variable name or column number of student identifier (ID) variable. ID variable must not contain any missing values.
}
  \item{wgt}{
Optional: Variable name or column number of weighting variable. If no weighting variable is specified,
all cases will be equally weighted.
}
  \item{type}{
Defines the replication method for cluster replicates which is to be applied. Depending on \code{type}, additional
arguments must be specified (e.g., \code{PSU} and/or \code{repInd} or \code{repWgt}).
}
  \item{PSU}{
Variable name or column number of variable indicating the primary sampling unit (PSU). When a jackknife procedure is applied,
the PSU is the jackknife zone variable. If \code{NULL}, no cluster structure is assumed and
standard errors are computed according to a random sample.
}
  \item{repInd}{
Variable name or column number of variable indicating replicate ID. In a jackknife procedure, this is the jackknife replicate
variable. If \code{NULL}, no cluster structure is assumed and standard errors are computed according to a random sample.
}
  \item{jkfac}{
Only applies if \code{engine = "BIFIEsurvey"}. Argument is passed to \code{\link[BIFIEsurvey]{BIFIE.data.jack}} and specifies
the factor for multiplying jackknife replicate weights.
}
  \item{repWgt}{
Normally, replicate weights are created by \code{repMean} directly from \code{PSU} and \code{repInd} variables. Alternatively,
if replicate weights are included in the data.frame, specify the variable names or column number in the \code{repWgt} argument.
}
\item{nest}{Optional: name or column number of the nesting variable. Only applies in nested multiple imputed data sets.
}
\item{imp}{Optional: name or column number of the imputation variable. Only applies in multiple imputed data sets.
}
\item{groups}{
Optional: vector of names or column numbers of one or more grouping variables.
}
  \item{group.splits}{
Optional: If groups are defined, \code{group.splits} optionally specifies whether analysis should be done also
in the whole group or overlying groups. See examples for more details.
}
  \item{group.differences.by}{
Optional: Specifies variable group differences should be computed for. The corresponding variable must be included in
the \code{groups} statement. Exception: choose 'wholePop' if you want to estimate each's group difference from the
overall sample mean. See examples for further details.
}
  \item{cross.differences}{
Either a list of vectors, specifying the pairs of levels for which cross-level differences should be computed.
Alternatively, if \code{TRUE}, cross-level differences for all pairs of levels are computed. If \code{FALSE}, no cross-level
differences are computed. (see example 2a, 3, and 4)
}
  \item{crossDiffSE}{
Method for standard error estimation for cross level differences, where groups are dependent.
\code{wec} uses weighted effect coding, \code{rep} uses replication methods (bootstrap or jackknife) to
estimate the standard error between the total mean and group-specific means. \code{old} does not account for dependent
groups and treat the groups as if they were independent from each other.
}
  \item{adjust}{
Variable name or column number of variable(s) for which adjusted means should be computed. Non-numeric variables (factors) will
be converted to 0/1 dichotomous variables.
}
  \item{useEffectLiteR}{
Logical: use the \code{lavaan}-wrapper \code{EffectLiteR} to compute adjusted means? Alternatively, adjusted means
are computed by applying a simple linear regression model in each group, using the variables in \code{adjust} as
independent variables. Afterwards, the coefficients are weighted with the (weighted) means of the independent
variables. Standard errors for this procedure are received using the delta method by applying an augmented
variance-covariance matrix which assumes zero covariances between independent variable means and regression
coefficients. We recommend to set \code{useEffectLiteR = TRUE} if no replication methods are applied. When
replication methods are used (jackknife-1, jackknife-2, BRR), we recommend to set \code{useEffectLiteR = FALSE},
because otherwise the estimation is very slow.
}
  \item{nBoot}{
Without replicates (i.e., for completely random samples), the \code{rep} method for standard error estimation for
cross level differences needs a bootstrap. \code{nBoot} therefore specifies the number of bootstrap samples.
This argument is only necessary, if \code{crossDiffSE = "rep"} \emph{and} none of the replicate methods
(JK1, JK2, or BRR) is applied. Otherwise, \code{nBoot} will be ignored.
}
  \item{group.delimiter}{
Character string which separates the group names in the output frame.
}
  \item{trend}{
Optional: name or column number of the trend variable which contains the measurement time of the survey.
Note: Levels of all grouping variables must be equal in all 'sub populations' partitioned by the discrete trend variable.
\code{repMean} computes differences for all pairwise contrasts defined by trend variable levels. or three measurement
occasions, i.e. 2010, 2015, and 2020, contrasts (i.e. trends) are computed for 2010 vs. 2015, 2010 vs. 2020, and
2015 vs. 2020.
}
  \item{linkErr}{
Optional: Either the name or column number of the linking error variable. If \code{NULL}, a linking error of 0 will be assumed in trend estimation.
Alternatively, linking errors may be
given as \code{data.frame} with following specifications: Two columns, named \code{trendLevel1} and \code{trendLevel2} which contain the
levels of the trend variable. The contrasts between both values indicates which trend is meant. For only two measurement occasions, i.e.
2010 and 2015, \code{trendLevel1} should be \code{2010}, and \code{trendLevel2} should be \code{2015}. For three measurement occasions,
i.e. 2010, 2015, and 2020, additional lines are necessary where \code{trendLevel1} should be \code{2010}, and \code{trendLevel2} should be
\code{2020}, to mark the contrast between 2010 and 2020, and further additional lines are necessary where \code{trendLevel1} should be
\code{2015}, and \code{trendLevel2} should be \code{2020}. The column \code{depVar} must include the name of the dependent variable. This
string must correspond to the name of the dependent variable in the data. The column \code{parameter} indicates the parameter the linking
error belongs to. Column \code{linkingError} includes the linking error value. Providing linking error in a data.frame is necessary for
more than two measurement occasions. See the example 3a for further details.
}
  \item{dependent}{
Variable name or column number of the dependent variable.
}
  \item{na.rm}{
Logical: Should cases with missing values be dropped?
}
  \item{doCheck}{
Logical: Check the data for consistency before analysis? If \code{TRUE} groups with insufficient data are excluded from
analysis to prevent subsequent functions from crashing.
}
  \item{engine}{
Which package should be used for estimation?
}
  \item{scale}{
scaling constant for variance, for details, see help page of \code{\link[survey]{svrepdesign}} from the \code{survey} package
}
  \item{rscales}{
scaling constant for variance, for details, see help page of \code{\link[survey]{svrepdesign}} from the \code{survey} package
}
  \item{mse}{
Logical: If \code{TRUE}, compute variances based on sum of squares around the point estimate, rather than the mean of the replicates.
See help page of \code{\link[survey]{svrepdesign}} from the \code{survey} package for further details.
}
  \item{rho}{
Shrinkage factor for weights in Fay's method. If \code{engine = "survey"}, argument is passed to the \code{rho} argument of
the \code{\link[survey]{svrepdesign}} function from the \code{survey} package. See the corresponding help page for further
details. If \code{engine = "BIFIEsurvey"}, argument is passed to the \code{fayfac} argument of the
\code{\link[BIFIEsurvey]{BIFIE.data.jack}} function from the \code{BIFIEsurvey} package. See the
corresponding help page for further details. For convenience, if \code{rho = NULL} (the default)
and \code{engine = "BIFIEsurvey"} and \code{type = "JK1"}, \code{\link[BIFIEsurvey]{BIFIE.data.jack}}
is called with \code{jktype="JK_GROUP"} and \code{fayfac = rho}, where \eqn{\rho = (N_{cluster} - 1) \times N_{cluster}^{-1}}
}
  \item{hetero}{
Logical: Assume heteroscedastic variance for weighted effect coding?
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
Software implementation used for estimating cross-level differences. Choices are either \code{"lavaan"} (required if
\code{stochasticGroupSites == "TRUE"}) or R function \code{\link[stats]{lm}}. \code{"lavaan"} is the default.
}
  \item{stochasticGroupSizes}{
Logical: Assume stochastic group sizes for using weighted effect coding in cross-level differences? Note: To date,
only \code{crossDiffSE.engine = "lavaan"} allows for stochastic group sizes. Stochastic group sizes are not yet
implemented for any replication method (jackknife, BRR).
}
  \item{verbose}{
Logical: Show analysis information on console?
}
  \item{progress}{
Logical: Show progress bar on console?
}
  \item{nCores}{
integer (default: NULL), number of cores to use for parallel processing, if \code{engine = "survey"}. If \code{NULL},
single core processing is used.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Function first creates replicate weights based on PSU and repInd variables (if defined) according to JK2 or BRR procedure
as implemented in WesVar. According to multiple imputed data sets, a workbook with several analyses is created.
The function afterwards serves as a wrapper for \code{\link[survey]{svymean}} called by \code{\link[survey]{svyby}} implemented in
the \sQuote{survey} package. The results of the several analyses are then pooled according to Rubin's rule.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of data frames in the long format. The output can be summarized using the \code{report} function.
The first element of the list is a list with either one (no trend analyses) or two (trend analyses)
data frames with at least six columns each. For each subpopulation denoted by the \code{groups}
statement, each parameter (i.e., mean, variance, or group differences) and each coefficient (i.e., the estimate 
and the corresponding standard error) the corresponding value is given.
\item{group}{Denotes the group an analysis belongs to. If no groups were specified and/or analysis for the 
whole sample were requested, the value of \sQuote{group} is \sQuote{wholeGroup}.}
\item{depVar}{Denotes the name of the dependent variable in the analysis. }
\item{modus}{Denotes the mode of the analysis. For example, if a JK2 analysis without sampling weights was 
conducted, \sQuote{modus} takes the value \sQuote{jk2.unweighted}. If a analysis without any replicates but with sampling
weights was conducted, \sQuote{modus} takes the value \sQuote{weighted}.}
\item{parameter}{Denotes the parameter of the regression model for which the corresponding value is given
further. Amongst others, the \sQuote{parameter} column takes the values \sQuote{mean}, \sQuote{sd}, \sQuote{var} and \sQuote{meanGroupDiff} if 
group differences were requested.}
\item{coefficient}{Denotes the coefficient for which the corresponding value is given further. Takes the 
values \sQuote{est} (estimate) and \sQuote{se} (standard error of the estimate).}
\item{value}{The value of the parameter estimate in the corresponding group.}
If groups were specified, further columns which are denoted by the group names are added to the data frame. 
}
\references{
te Grotenhuis, M., Pelzer, B., Eisinga, R., Nieuwenhuis, R., Schmidt-Catran, A., & Konig, R. (2017).
When size matters: advantages of weighted effect coding in observational studies. \emph{International Journal of Public Health.} \bold{62}, 163--167.

Sachse, K. A. & Haag, N. (2017). Standard errors for national trends in international
large-scale assessments in the case of cross-national differential item functioning. \emph{Applied
Measurement in Education, 30}, (2), 102-116. http://dx.doi.org/10.1080/08957347.2017.1283315

Weirich, S., Hecht, M., Becker, B. et al. Comparing group means with the total mean in random samples,
surveys, and large-scale assessments: A tutorial and software illustration. Behav Res (2021).
https://doi.org/10.3758/s13428-021-01553-1
}
\examples{
data(lsa)

### Example 1: only means, SD and variances for each country
### We only consider domain 'reading'
rd     <- lsa[which(lsa[,"domain"] == "reading"),]

### We only consider the first "nest".
rdN1   <- rd[which(rd[,"nest"] == 1),]

### First, we only consider year 2010
rdN1y10<- rdN1[which(rdN1[,"year"] == 2010),]

\donttest{
### mean estimation
means1 <- repMean(datL = rdN1y10, ID="idstud", wgt="wgt", type = "JK2", PSU = "jkzone",
          repInd = "jkrep", imp="imp", groups = "country", dependent = "score",
          na.rm=FALSE, doCheck=TRUE, engine = "BIFIEsurvey")
### reporting function: the function does not know which content domain is being considered,
### so the user may add new columns in the output using the 'add' argument
res1   <- report(means1, add = list(domain = "reading"))

### Example 1a: Additionally to example 1, we decide to estimate whether
### each country's mean differ significantly from the overall mean as well
### as from the individual means of the other contries
means1a<- repMean(datL = rdN1y10, ID="idstud", wgt="wgt", type = "JK2", PSU = "jkzone",
          repInd = "jkrep", imp="imp", groups = "country", group.splits = 0:1,
          group.differences.by = "country", cross.differences = TRUE,
          dependent = "score", na.rm=FALSE, doCheck=TRUE, hetero=FALSE)
res1a  <- report(means1a, add = list(domain = "reading"))

### See that the means of all countries significantly differ from the overall mean.
print(res1a[intersect(which(res1a[,"comparison"] == "crossDiff"),
      which(res1a[,"parameter"] == "mean")),], digits = 3)

### Example 2: Sex differences by country. Assume equally weighted cases by omitting
### 'wgt' argument.
means2 <- repMean(datL = rdN1y10, ID="idstud", type = "JK2", PSU = "jkzone",
          repInd = "jkrep", imp="imp", groups = c("country", "sex"), group.splits = 0:2,
          group.differences.by="sex", dependent = "score", na.rm=FALSE, doCheck=TRUE,
          cross.differences =TRUE, crossDiffSE.engine= "lm")
res2   <- report(means2,add = list(domain = "reading"))

### Example 2a: Additionally to example 2, we decide to estimate whether
### each country's mean differ significantly from the overall mean. (Note: by default,
### such cross level differences are estimated using 'weighted effect coding'. Use the
### 'crossDiffSE' argument to choose alternative methods.) Moreover, we estimate whether
### each country's sex difference differ significantly from the sex difference in the
### whole population.
means2a<- repMean(datL = rdN1y10, ID="idstud", wgt="wgt", type = "JK2", PSU = "jkzone",
          repInd = "jkrep", imp="imp", groups = c("country", "sex"), group.splits = 0:2,
          group.differences.by="sex", cross.differences = list(c(0,1), c(0,2)),
          dependent = "score", na.rm=FALSE, doCheck=TRUE,
          crossDiffSE.engine= "lm", clusters = "idclass")
res2a  <- report(means2a,add = list(domain = "reading"))

### Third example: like example 2a, but using nested imputations of dependent variable,
### and additionally estimating trend: use 'rd' instead of 'rdN1y10'
### assume equally weighted cases by omitting 'wgt' argument
### ignoring jackknife by omitting 'type', 'PSU' and 'repInd' argument
means3T<- repMean(datL = rd, ID="idstud", imp="imp", nest="nest",
          groups = c("country", "sex"), group.splits = 0:2, group.differences.by="sex",
          cross.differences = list(c(0,1), c(0,2)), dependent = "score", na.rm=FALSE,
          doCheck=TRUE, trend = "year", linkErr = "leScore",
          crossDiffSE = "wec", crossDiffSE.engine= "lavaan")
res3T  <- report(means3T, add = list(domain = "reading"))

### Example 3a: like example 3, but providing linking errors in an additional data.frame
### This is optional for two measurement occasions but mandatory if the analysis contains
### more than two measurement occasions
linkErr<- data.frame ( trendLevel1 = 2010, trendLevel2 = 2015,  depVar = "score",
          parameter = "mean", unique(lsa[,c("domain", "leScore")]),
          stringsAsFactors = FALSE)
colnames(linkErr) <- car::recode(colnames(linkErr), "'leScore'='linkingError'")
### note that the linking errors for the specified domain have to be chosen via
### subsetting
means3a<- repMean(datL = rd, ID="idstud", imp="imp", nest="nest",
          groups = c("country", "sex"),
          group.splits = 0:2, group.differences.by="sex",
          cross.differences = list(c(0,1), c(0,2)),
          dependent = "score", na.rm=FALSE, doCheck=TRUE, trend = "year",
          linkErr = linkErr[which(linkErr[,"domain"] == "reading"),],
          crossDiffSE = "wec", crossDiffSE.engine= "lavaan")
res3a  <- report(means3a, add = list(domain = "reading"))

### Fourth example: using a loop do analyse 'reading' and 'listening' comprehension
### in one function call. Again with group and cross differences and trends, and
### trend differences
### we use weights but omit jackknife analysis by omitting 'type', 'PSU' and 'repInd'
### argument
means4T<- by ( data = lsa, INDICES = lsa[,"domain"], FUN = function (sub.dat) {
          repMean(datL = sub.dat, ID="idstud", wgt="wgt", imp="imp", nest="nest",
                 groups = c("country", "sex"), group.splits = 0:2,
                 group.differences.by="sex",
                 cross.differences = list(c(0,1), c(0,2)), dependent = "score",
                 na.rm=FALSE, doCheck=TRUE,
                 trend = "year", linkErr = "leScore", crossDiffSE.engine= "lm") })
ret4T  <- do.call("rbind", lapply(names(means4T), FUN = function ( domain ) {
          report(means4T[[domain]], add = list(domain = domain))}))
          
### Fifth example: compute adjusted means, also with trend estimation
### Note: all covariates must be numeric or 0/1 dichotomous
rdN1[,"mignum"] <- as.numeric(rdN1[,"mig"])
rdN1[,"sexnum"] <- car::recode(rdN1[,"sex"], "'male'=0; 'female'=1", as.numeric=TRUE,
                   as.factor=FALSE)
means5 <- repMean(datL = rdN1, ID="idstud", wgt="wgt", type = "JK2", PSU = "jkzone",
          repInd = "jkrep", imp="imp", groups = "country",
          adjust = c("sexnum", "ses", "mignum"), useEffectLiteR = FALSE,
          dependent = "score", na.rm=FALSE, doCheck=TRUE, trend = "year",
          linkErr = "leScore")
res5   <- report(means5, add = list(domain = "reading"))

\dontrun{
############################################################################################
#    Example 6: R code for running the PISA 2015 science example to compare group means    #
#                    with the total mean using weighted effect coding                      #
############################################################################################

# Warning: large PISA data set requires at least 16 GB free working memory (RAM):

### define necessary directories (note: writing permissions required)
folder <- tempdir()

### download PISA 2015 zipped student questionnaire data (420 MB) to a folder with
### writing permissions
download.file(url = "https://webfs.oecd.org/pisa/PUF_SPSS_COMBINED_CMB_STU_QQQ.zip",
         destfile = file.path(folder, "pisa2015.zip"))

### unzip PISA 2015 student questionnaire data (1.5 GB) to temporary folder
zip::unzip(zipfile = file.path(folder, "pisa2015.zip"), files= "CY6_MS_CMB_STU_QQQ.sav",
     exdir=folder)

### read data
pisa <- foreign::read.spss(file.path (folder, "CY6_MS_CMB_STU_QQQ.sav"),
        to.data.frame=TRUE, use.value.labels = FALSE, use.missings = TRUE)

# dependent variables
measure.vars <- paste0("PV", 1:10, "SCIE")

### choose desired variables and reshape into the long format
#              'CNTSTUID' = individual student identifier
#                   'CNT' = country identifier
#                 'SENWT' = senate weight (assume a population of 5000 in each country)
#              'W_FSTUWT' = final student weight
#                  'OECD' = dummy variable indicating which country is part of the OECD
#   'W_FSTURWT' (1 to 80) = balanced repeated replicate weights
# 'PV1SCIE' to 'PV10SCIE' = 10 plausible values of (latent) science performance
pisaLong <- reshape2::melt(pisa, id.vars = c("CNTSTUID", "CNT", "SENWT", "W_FSTUWT",
            "OECD", paste0("W_FSTURWT", 1:80)),
            measure.vars = measure.vars, value.name = "value", variable.name="imp",
            na.rm=TRUE)

### choose OECD countries
oecd <- pisaLong[which(pisaLong[,"OECD"] == 1),]

### analyze data
### analysis takes approximately 30 minutes on an Intel i5-6500 machine with 32 GB RAM
means   <- repMean( datL = oecd,         # data.frame in the long format
    ID                   = "CNTSTUID",   # student identifier
    dependent            = "value",      # the dependent variable in the data
    groups               = "CNT",        # the grouping variable
    wgt                  = "SENWT",      # (optional) weighting variable. We use senate
                                         # weights (assume a population of 5000 in each
                                         # country)
    type                 = "Fay",        # type of replication method. Corresponding to
                                         # the PISA sampling method, we use "Fay"
    rho                  = 0.5,          # shrinkage factor for weights in Fay's method
    scale                = NULL,         # scaling constant for variance, set to NULL
                                         # according to PISA's sampling method
    rscales              = NULL,         # scaling constant for variance, set to NULL
                                         # according to PISA's sampling method
    repWgt               = paste0("W_FSTURWT", 1:80), # the replicate weights,
                                                      # provided by the OECD
    imp                  = "imp",        # the imputation variable
    mse                  = FALSE,        # if TRUE, compute variances based on sum of
                                         # squares around the point estimate, rather
                                         # than the mean of the replicates.
    group.splits         = 0:1,          # defining the 'levels' for which means should
                                         # be computed. 0:1 implies that means for the
                                         # whole sample (level 0) as well as for groups
                                         # (level 1) are computed
    cross.differences    = TRUE,         # defines whether (and which) cross level mean
                                         # differences should be computed. TRUE means
                                         # that all cross level mean differences are
                                         # computed
    crossDiffSE          = "wec",        # method for standard errors of mean
                                         # differences
    crossDiffSE.engine   = "lm",         # software implementation for standard
                                         # errors of mean differences
    hetero               = TRUE,         # assume heteroscedastic group variances
    stochasticGroupSizes = FALSE         # assume fixed group sizes
                   )

### call a reporting function to generate user-friendly output
results <- report(means, exclude = c("Ncases", "NcasesValid", "var", "sd"))
}}}


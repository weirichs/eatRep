\name{repTable}
\alias{repTable}
\alias{jk2.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{JK1, JK2 and BRR for frequency tables and trend estimation.}
\description{Compute frequency tables for categorical variables (e.g. factors: dichotomous or polytomous) in complex
cluster designs. Estimation of standard errors optionally takes the clustered structure and multiple imputed
variables into account. To date, Jackknife-1 (JK1), Jackknife-2 (JK2) and Balanced repeated replicate (BRR) methods are 
implemented to account for clustered designs. Procedures of Rubin (1987) and Rubin (2003) are implemented to account for 
multiple imputed data and nested imputed data, if necessary. Conceptually, the function combines replication and imputation 
methods. Technically, this is a wrapper for the \code{\link[survey]{svymean}} function of the \code{survey} package.
}
\usage{
repTable(datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"), PSU = NULL,
          repInd = NULL, jkfac=NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL,
          group.splits = length(groups), group.differences.by = NULL,
          cross.differences = FALSE, crossDiffSE = c("wec", "rep","old"),
          nBoot = 100, chiSquare = FALSE, correct = TRUE, group.delimiter = "_",
          trend = NULL, linkErr = NULL, dependent, separate.missing.indicator = FALSE,
          na.rm=FALSE, expected.values = NULL, doCheck = TRUE, forceTable = FALSE,
          engine = c("survey", "BIFIEsurvey"), scale = 1, rscales = 1, mse=TRUE,
          rho=NULL, verbose = TRUE, progress = TRUE )}
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
Normally, replicate weights are created by \code{repTable} directly from \code{PSU} and \code{repInd} variables. Alternatively,
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
Optional: Specifies one grouping variable for which a chi-square test should be applied.
The corresponding variable must be included in the \code{groups} statement. If specified, the
distribution of the dependent variable is compared between the groups. See examples for
further details.
}
  \item{cross.differences}{
Either a list of vectors, specifying the pairs of levels for which cross-level differences should be computed.
Alternatively, if TRUE, cross-level differences for all pairs of levels are computed. If FALSE, no cross-level
differences are computed. (see examples 2a, 3, and 4 in the help file of the \code{repMean} function)
}
  \item{crossDiffSE}{
Method for standard error estimation for cross level differences, where groups are dependent.
\code{wec} uses weighted effect coding, \code{rep} uses replication methods (bootstrap or jackknife) to
estimate the standard error between the total mean and group-specific means. \code{old} does not account for dependent
groups and treat the groups as if they were independent from each other.
}
  \item{nBoot}{
Without replicates (i.e., for completely random samples), the \code{rep} method for standard error estimation for
cross level differences needs a bootstrap. \code{nBoot} therefore specifies the number of bootstrap samples.
This argument is only necessary, if \code{crossDiffSE = "rep"} \emph{and} none of the replicate methods
(JK1, JK2, or BRR) is applied. Otherwise, \code{nBoot} will be ignored.
}
  \item{chiSquare}{
Logical. Applies only if \code{group.differences.by} was specified. Defines whether
group differences should be represented in a chi square test or in (mean) differences of each
group's relative frequency. Note: To date, chi square test is not available for \code{engine = "BIFIEsurvey"}.
}
  \item{correct}{
Logical. Applies only if 'group.differences.by' is requested without cluster replicates. A logical indicating whether to apply
continuity correction when computing the test statistic for 2 by 2 tables. See help page of 'chisq.test' for further details.
}
  \item{group.delimiter}{
Character string which separates the group names in the output frame.
}
  \item{trend}{
Optional: name or column number of the trend variable. Note: Trend variable must have exact two levels. Levels for
grouping variables must be equal in both 'sub populations' partitioned by the trend variable.
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
more than two measurement occasions. See the fourth example below for further details.
}
  \item{dependent}{
Variable name or column number of the dependent variable.
}
  \item{separate.missing.indicator}{
Logical. Should frequencies of missings in dependent variable be integrated? Note: That is only useful if missing occur as \code{NA}. If the dependent variable
is coded as character, for example \code{'male', 'female', 'missing'}, separate missing indicator is not necessary.
}
  \item{na.rm}{
Logical: Should cases with missing values be dropped?
}
  \item{expected.values}{
Optional. A vector of values expected in dependent variable. Recommend to left this argument empty.
}
  \item{doCheck}{
Logical: Check the data for consistency before analysis? If \code{TRUE} groups with insufficient data are excluded from
analysis to prevent subsequent functions from crashing.
}
  \item{forceTable}{
Logical: Function decides internally whether the table or the mean function of \code{survey} is called.
If the mean function is called, the polytomous dependent variable is converted to dichotomous indicator
variables. If mean is called, group differences for each category of the polytomous dependent variable
can be computed. If table is called, a chi square statistic may be computed. The argument allows to
force the function either to call mean or table.
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
  \item{verbose}{
Logical: Show analysis information on console?
}
  \item{progress}{
Logical: Show progress bar on console?
}
}
\details{
Function first creates replicate weights based on PSU and repInd variables according to JK2 procedure
implemented in WesVar. According to multiple imputed data sets, a workbook with several analyses is created.
The function afterwards serves as a wrapper for \code{\link[survey]{svymean}} called by \code{\link[survey]{svyby}} implemented in the \code{survey} package.
Relative frequencies of the categories of the dependent variable are computed by the means of the dichotomous indicators
(e.g. dummy variables) of each category. The results of the several analyses are then pooled according to Rubin's rule,
which is adapted for nested imputations if the \code{dependent} argument implies a nested structure.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A list of data frames in the long format. The output can be summarized using the \code{report} function.
The first element of the list is a list with either one (no trend analyses) or two (trend analyses)
data frames with at least six columns each. For each subpopulation denoted by the \code{groups} statement, each
dependent variable, each parameter (i.e., the values of the corresponding categories of the dependent variable)
and each coefficient (i.e., the estimate and the corresponding standard error) the corresponding value is given.
\item{group}{Denotes the group an analysis belongs to. If no groups were specified and/or analysis for the 
whole sample were requested, the value of \sQuote{group} is \sQuote{wholeGroup}.}
\item{depVar}{Denotes the name of the dependent variable in the analysis. }
\item{modus}{Denotes the mode of the analysis. For example, if a JK2 analysis without sampling weights was 
conducted, \sQuote{modus} takes the value \sQuote{jk2.unweighted}. If a analysis without any replicates but with sampling
weights was conducted, \sQuote{modus} takes the value \sQuote{weighted}.}
\item{parameter}{Denotes the parameter of the regression model for which the corresponding value is given
further. For frequency tables, this is the value of the category of the dependent variable which relative
frequency is given further.}
\item{coefficient}{Denotes the coefficient for which the corresponding value is given further. Takes the 
values \sQuote{est} (estimate) and \sQuote{se} (standard error of the estimate).}
\item{value}{The value of the parameter, i.e. the relative frequency or its standard error.}
If groups were specified, further columns which are denoted by the group names are added to the data frame. 
}
\references{
  Rubin, D.B. (2003): Nested multiple imputation of NMES
  via partially incompatible MCMC.
  \emph{Statistica Neerlandica} \bold{57, 1}, 3--18.
}
\examples{
data(lsa)

### Example 1: only means, SD and variances for each country
### subsetting: We only consider domain 'reading'
rd     <- lsa[which(lsa[,"domain"] == "reading"),]

### We only consider the first "nest".
rdN1   <- rd[which(rd[,"nest"] == 1),]

### First, we only consider year 2010
rdN1y10<- rdN1[which(rdN1[,"year"] == 2010),]

### First example: Computes frequencies of polytomous competence levels (1, 2, 3, 4, 5)
### conditionally on country, using a chi-square test to decide whether the distribution
### varies between countries (it's an overall test, i.e. with three groups, df1=8).
freq.tab1 <- repTable(datL = rdN1y10, ID = "idstud", wgt = "wgt", imp="imp",
             type = "JK2", PSU = "jkzone", repInd = "jkrep", groups = "country",
             group.differences.by = "country", dependent = "comp", chiSquare = TRUE)
res1      <- report(freq.tab1, add = list ( domain = "reading" ))

\donttest{
### Second example: Computes frequencies of polytomous competence levels (1, 2, 3, 4, 5)
### conditionally on country. Now we test whether the frequency of each single category
### differs between pairs of countries (it's not an overall test ... repTable now
### calls repMean internally, using dummy variables
freq.tab2 <- repTable(datL = rdN1y10, ID = "idstud", wgt = "wgt", imp="imp",
             type = "JK2", PSU = "jkzone", repInd = "jkrep", groups = "country",
             group.differences.by = "country", dependent = "comp", chiSquare = FALSE)
res2      <- report(freq.tab2, add = list ( domain = "reading" ))

### Third example: trend estimation and nested imputation and 'by' loop
### (to date, only crossDiffSE = "old" works)
freq.tab3 <- by ( data = lsa, INDICES = lsa[,"domain"], FUN = function (subdat) {
             repTable(datL = subdat, ID = "idstud", wgt = "wgt", imp="imp",
                 nest = "nest", type = "JK2", PSU = "jkzone", repInd = "jkrep",
                 groups = "country", group.differences.by = "country",
                 group.splits = 0:1, cross.differences = TRUE, crossDiffSE = "old",
                 dependent = "comp", chiSquare = FALSE, trend = "year",
                 linkErr = "leComp") })
res3      <- do.call("rbind", lapply(names(freq.tab3), FUN = function (domain) {
             report(freq.tab3[[domain]], trendDiffs = TRUE,
                    add = list ( domain = domain ))
             }))
             
### Fourth example: similar to example 3. trend estimation using a linking
### error data.frame
linkErrs  <- data.frame ( trendLevel1 = 2010, trendLevel2 = 2015,  depVar = "comp",
             unique(lsa[,c("domain", "comp", "leComp")]), stringsAsFactors = FALSE)
colnames(linkErrs) <- car::recode(colnames(linkErrs),
                      "'comp'='parameter'; 'leComp'='linkingError'")
freq.tab4 <- by ( data = lsa, INDICES = lsa[,"domain"], FUN = function (subdat) {
             repTable(datL = subdat, ID = "idstud", wgt = "wgt", type="none",
                 imp="imp", nest = "nest", groups = "country",
                 group.differences.by = "country", group.splits = 0:1,
                 cross.differences = FALSE, dependent = "comp", chiSquare = FALSE,
                 trend = "year",
                 linkErr = linkErrs[which(linkErrs[,"domain"] == subdat[1,"domain"]),])
             })
res4      <- do.call("rbind", lapply(names(freq.tab4), FUN = function (domain) {
             report(freq.tab4[[domain]], trendDiffs = TRUE,
                    add = list ( domain = domain ))
             }))

### Fifth example: minimal example for three measurement occasions
### borrow data from the eatGADS package
trenddat1 <- system.file("extdata", "trend_gads_2010.db", package = "eatGADS")
trenddat2 <- system.file("extdata", "trend_gads_2015.db", package = "eatGADS")
trenddat3 <- system.file("extdata", "trend_gads_2020.db", package = "eatGADS")
trenddat  <- eatGADS::getTrendGADS(filePaths = c(trenddat1, trenddat2, trenddat3),
             years = c(2010, 2015, 2020), fast=FALSE)
dat       <- eatGADS::extractData(trenddat)
### use template linking Error Object
load(system.file("extdata", "linking_error.rda", package = "eatRep"))
### check consistency of data and linking error object
check1 <- checkLEs(c(trenddat1, trenddat2, trenddat3), lErr)
### Analysis for reading comprehension
freq.tab5 <- repTable(datL = dat[which(dat[,"dimension"] == "reading"),],
             ID = "idstud", type="none", imp="imp", dependent = "traitLevel",
             chiSquare = FALSE, trend = "year",
             linkErr = lErr[which(lErr[,"domain"] == "reading"),])
res5      <- report(freq.tab5, trendDiffs = TRUE, add = list ( domain = "reading" ))
}}

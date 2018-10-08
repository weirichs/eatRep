\name{jk2.mean}
\alias{jk2.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{JK1, JK2 and BRR for mean estimates.}
\description{Compute totals, means, variances and standard deviations with standard errors
for complex cluster designs with multiple imputed variables (e.g. plausible values) based on Jackknife (JK1, JK2)
or Balanced Repeated Replicates (BRR) procedure. Conceptually, the function combines replication methods and 
methods for multiple imputed data. Nested imputations of the dependent variable(s) are supported as well.
Technically, this is a wrapper for the \code{svymean()} and \code{svyvar()} functions of the \code{'survey'} package.}
\usage{
jk2.mean (datL, ID, wgt = NULL, type = c("JK1", "JK2", "BRR"), PSU = NULL, repInd = NULL, 
          repWgt = NULL, nest=NULL, imp=NULL, groups = NULL, group.splits = length(groups),
            group.differences.by = NULL, cross.differences = FALSE, group.delimiter = "_",
            trend = NULL, linkErr = NULL, dependent, na.rm = FALSE, doCheck = TRUE)}
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
Defines the replication method for cluster replicates which is to be applied. Without cluster replicates (i.e., if \code{PSU} 
and/or \code{repInd} is NULL, \code{type} will be ignored.
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
Normally, replicate weights are created by \code{jk2.mean} directly from \code{PSU} and \code{repInd} variables. Alternatively, 
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
  \item{group.differences.by}{
%%     ~~Describe \code{sig.dif.bound} here~~
Optional: Specifies variable group differences should be computed for. The corresponding variable must be included in 
the \code{groups} statement. Exception: choose 'wholePop' if you want to estimate each's group difference from the 
overall sample mean. See examples for further details.
}
  \item{cross.differences}{
%%     ~~Describe \code{sig.dif.bound} here~~
Either a list of vectors, specifying the pairs of levels for which cross-level differences should be computed.
Alternatively, if TRUE, cross-level differences for all pairs of levels are computed. If FALSE, no cross-level
differences are computed. (see example 2a, 3, and 4)
}
  \item{group.delimiter}{
%%     ~~Describe \code{sig.dif.bound} here~~
Character string which separates the group names in the output frame.
}
  \item{trend}{
%%     ~~Describe \code{sig.dif.bound} here~~
Optional: name or column number of the trend variable. Note: Trend variable must have exact two levels. Levels for 
grouping variables must be equal in both 'sub populations' partitioned by the trend variable. 
}
  \item{linkErr}{
%%     ~~Describe \code{sig.dif.bound} here~~
Optional: name or column number of the trend variable. If 'NULL', a linking error of 0 will be assumed in trend estimation. 
}
  \item{dependent}{
%%     ~~Describe \code{sig.dif.bound} here~~
Variable name or column number of the dependent variable.
}
  \item{na.rm}{
Logical: Should cases with missing values be dropped?
}
  \item{doCheck}{
Logical: Check the data for consistency before analysis? If \code{TRUE} groups with insufficient data are excluded from
analysis to prevent subsequent functions from crashing.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Function first creates replicate weights based on PSU and repInd variables (if defined) according to JK2 or BRR procedure
as implemented in WesVar. According to multiple imputed data sets, a workbook with several analyses is created.
The function afterwards serves as a wrapper for \code{svymean()} called by \code{svyby()} implemented in
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
%% ~put references to the literature/web site here ~
}
\author{
Sebastian Weirich
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

### mean estimation
means1 <- jk2.mean(datL = rdN1y10, ID="idstud", wgt="wgt", type = "JK2", PSU = "jkzone", repInd = "jkrep",
          imp="imp", groups = "country", dependent = "score", na.rm=FALSE, doCheck=TRUE)
### reporting function: the function does not know which content domain is being considered,
### so it is possible to add new columns in the output using the 'add' argument
res1   <- report(means1, add = list(domain = "reading"))

### Example 1a: Additionally to example 1, we decide to estimate whether 
### each country's mean differ significantly from the overall mean as well
### as from the individual means of the other contries
means1a<- jk2.mean(datL = rdN1y10, ID="idstud", wgt="wgt", type = "JK2", PSU = "jkzone", repInd = "jkrep",
          imp="imp", groups = "country", group.splits = 0:1, group.differences.by = "country",
          cross.differences = TRUE, dependent = "score", na.rm=FALSE, doCheck=TRUE)
res1a  <- report(means1a, add = list(domain = "reading"))

### See that only the mean of 'LandB' significantly differs from the overall mean.
print(res1a[intersect(grep("wholeGroup.vs.", res1a[,"group"]), which(res1a[,"parameter"] == "mean")), ], digits = 3)

### Example 2: Sex differences by country. Assume equally weighted cases by omitting
### 'wgt' argument.
means2 <- jk2.mean(datL = rdN1y10, ID="idstud", type = "JK2", PSU = "jkzone", repInd = "jkrep",
          imp="imp", groups = c("country", "sex"), group.splits = 0:2, group.differences.by="sex",
          dependent = "score", na.rm=FALSE, doCheck=TRUE)
res2   <- report(means2,add = list(domain = "reading"))

### Example 2a: Additionally to example 2, we decide to estimate whether
### each country's mean differ significantly from the overall mean. Moreover, we estimate
### whether each country's sex difference differ significantly from the sex difference in
### the whole population. 
means2a<- jk2.mean(datL = rdN1y10, ID="idstud", type = "JK2", PSU = "jkzone", repInd = "jkrep",
          imp="imp", groups = c("country", "sex"), group.splits = 0:2, group.differences.by="sex",
          cross.differences = list(c(0,1), c(0,2)), dependent = "score", na.rm=FALSE, doCheck=TRUE)
res2a  <- report(means2a,add = list(domain = "reading"))

### Third example: like example 2a, but using nested imputations of dependent variable,
### and additionally estimating trend: use 'rd' indtead of 'rdN1y10'
means3T<- jk2.mean(datL = rd, ID="idstud", type = "JK2", PSU = "jkzone", repInd = "jkrep",
          imp="imp", nest="nest", groups = c("country", "sex"), group.splits = 0:2, group.differences.by="sex",
          cross.differences = list(c(0,1), c(0,2)), dependent = "score", na.rm=FALSE, doCheck=TRUE,
          trend = "year", linkErr = "leScore")
res3T  <- report(means3T, add = list(domain = "reading"))
          
### Fourth example: using a loop do analyse 'reading' and 'listening' comprehension
### in one function call. Again with group and cross differences and trends, and
### trend differences
means4T<- by ( data = lsa, INDICES = lsa[,"domain"], FUN = function (sub.dat) {
          jk2.mean(datL = sub.dat, ID="idstud", type = "JK2", PSU = "jkzone", repInd = "jkrep",
                 imp="imp", nest="nest", groups = c("country", "sex"), group.splits = 0:2, group.differences.by="sex",
                 cross.differences = list(c(0,1), c(0,2)), dependent = "score", na.rm=FALSE, doCheck=TRUE,
                 trend = "year", linkErr = "leScore") })
ret4T  <- do.call("rbind", lapply(names(means4T), FUN = function ( domain ) {
          report(means4T[[domain]], trendDiffs = TRUE, add = list(domain = domain))}))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

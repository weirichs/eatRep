\name{jk2.quantile}
\alias{jk2.quantile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{JK1, JK2 and BRR method for quantiles and trend estimation.}
\description{Compute quantiles with standard errors for complex cluster designs with multiple imputed variables
(e.g. plausible values) based on Jackknife (JK1, JK2) or balanced repeated replicates (BRR) procedure. Conceptually, 
the function combines replication methods and methods for multiple imputed data. Technically, this is a wrapper for
the \code{svyquantile()} function of the survey package.
}
\usage{
jk2.quantile(datL, ID, wgt = NULL, type = c("JK1", "JK2", "BRR"), PSU = NULL, repInd = NULL, 
            repWgt = NULL, nest=NULL, imp=NULL, groups = NULL, group.splits = length(groups),
            cross.differences = FALSE, group.delimiter = "_", trend = NULL, linkErr = NULL, dependent,
            probs = seq(0, 1, 0.25),  na.rm = FALSE, nBoot = NULL, 
            bootMethod = c("wSampling","wQuantiles") , doCheck = TRUE)}
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
Normally, replicate weights are created by \code{jk2.quantile} directly from \code{PSU} and \code{repInd} variables. Alternatively, 
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
  \item{cross.differences}{
%%     ~~Describe \code{sig.dif.bound} here~~
Either a list of vectors, specifying the pairs of levels for which cross-level differences should be computed.
Alternatively, if TRUE, cross-level differences for all pairs of levels are computed. If FALSE, no cross-level
differences are computed. (see examples 2a, 3, and 4 in the help file of jk2.mean)
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
Alternatively, the linking error may be given as a single scalar value (i.e. 'linkErr = 1.225'). 
}
  \item{dependent}{
%%     ~~Describe \code{sig.dif.bound} here~~
Variable name or column number of the dependent variable.
}
  \item{probs}{
%%     ~~Describe \code{sig.dif.bound} here~~
Numeric vector with probabilities for which to compute quantiles.
}
  \item{na.rm}{
%%     ~~Describe \code{sig.dif.bound} here~~
Logical: Should cases with missing values be dropped?
}
  \item{nBoot}{
%%     ~~Describe \code{sig.dif.bound} here~~
Optional: Without replicates, standard error cannot be computed in a weighted sample. Alternatively, standard errors may
be computed using the \code{boot} package. \code{nBoot} therefore specifies the number of bootstrap samples. If not specified,
no standard errors will be given. In analyses containing replicates or samples without specifying person weights,
\code{nBoot} will be ignored.
}
  \item{bootMethod}{
%%     ~~Describe \code{sig.dif.bound} here~~
Optional: If standard error are computed in a bootstrap, two possible methods may be applied.
\code{'wSampling'} requests the function to draw \code{nBoot} weighted bootstrap samples for which unweighted quantiles
are computed. \code{'wQuantiles'} requests the function to draw \code{nBoot} unweighted bootstrap samples for which
weighted quantiles are computed.
}
  \item{doCheck}{
Logical: Check the data for consistency before analysis? If \code{TRUE} groups with insufficient data are excluded from
analysis to prevent subsequent functions from crashing.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Function first creates replicate weights based on PSU and repInd variables according to JK2 or BRR procedure
implemented in WesVar. According to multiple imputed data sets, a workbook with several analyses is created.
The function afterwards serves as a wrapper for \code{svyquantile()} called by \code{svyby()} implemented in the 'survey' package.
The results of the several analyses are then pooled according to Rubins rule, which is adapted for nested imputations if the \code{dependent} argument implies a nested structure.
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
%% ~put references to the literature/web site here ~
}
\author{
Sebastian Weirich
}
\examples{
### Example 1: only means, SD and variances for each country
### We only consider domain 'reading'
rd     <- lsa[which(lsa[,"domain"] == "reading"),]

### We only consider the first "nest".
rdN1   <- rd[which(rd[,"nest"] == 1),]

### First, we only consider year 2010
rdN1y10<- rdN1[which(rdN1[,"year"] == 2010),]

### First example: Computes percentile in a nested data structure for reading 
### scores conditionally on country and for the whole group 
perzent   <- jk2.quantile(datL = rd, ID = "idstud", wgt = "wgt", type = "JK2",
             PSU = "jkzone", repInd = "jkrep", imp = "imp", nest="nest",
             groups = "country", group.splits = c(0:1), dependent = "score", 
             probs = seq(0.1,0.9,0.2) )
res       <- report(perzent, add = list(domain = "reading"))

### Second example: Computes percentile for reading scores conditionally on country,
### use 100 bootstrap samples, assume no nested structure 
perzent   <- jk2.quantile(datL = rdN1y10, ID = "idstud", wgt = "wgt",
             imp = "imp", groups = "country", dependent = "score",
             probs = seq(0.1,0.9,0.2), nBoot = 100 )
res       <- report(perzent, add = list(domain = "reading"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

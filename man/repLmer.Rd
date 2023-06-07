\name{repLmer}
\alias{repLmer}
\title{Replication methods (JK1 and JK2) for multilevel linear regression models and trend estimation.}
\description{Compute multilevel linear models for complex cluster designs with multiple imputed variables based
on the Jackknife (JK1, JK2) procedure. Conceptually, the function combines replication
methods and methods for multiple imputed data. Technically, this is a wrapper for the \code{\link[BIFIEsurvey]{BIFIE.twolevelreg}} function
of the \code{BIFIEsurvey} package. \code{repLmer} only adds functionality for trend estimation. Please note
that the function is not suitable for logistic logit/probit models.
}
\usage{
repLmer(datL, ID, wgt = NULL, L1wgt=NULL, L2wgt=NULL, type = c("JK2", "JK1"),
            PSU = NULL, repInd = NULL, jkfac = NULL, rho = NULL, imp=NULL,
            group = NULL, trend = NULL, dependent, formula.fixed, formula.random,
            doCheck = TRUE, na.rm = FALSE, clusters, verbose = TRUE)}
\arguments{
  \item{datL}{
Data frame in the long format (i.e. each line represents one ID unit in one imputation of one nest) containing all
variables for analysis.
}
  \item{ID}{
Variable name or column number of student identifier (ID) variable. ID variable must not contain any missing values.
}
  \item{wgt}{
Optional: Variable name or column number of case weighting variable. If no weighting variable is specified,
all cases will be equally weighted.
}
  \item{L1wgt}{
Name of Level 1 weight variable. This is optional. If it is not provided, \code{L1wgt} is calculated from the
total weight (i.e., \code{wgt}) and \code{L2wgt}.
}
  \item{L2wgt}{
Name of Level 2 weight variable
}
  \item{type}{
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
  \item{jkfac}{
Argument is passed to \code{\link[BIFIEsurvey]{BIFIE.data.jack}} and specifies the factor for multiplying jackknife replicate
weights.
}
  \item{rho}{
Fay factor for statistical inference. The argument is passed to the \code{fayfac} argument of the
\code{\link[BIFIEsurvey]{BIFIE.data.jack}} function from the \code{BIFIEsurvey} package. See the
corresponding help page for further details. For convenience, if \code{rho = NULL} (the default)
and \code{type = "JK1"}, \code{\link[BIFIEsurvey]{BIFIE.data.jack}} is called with
\code{jktype="JK_GROUP"} and \code{fayfac = rho}, where \eqn{\rho = (N_{cluster} - 1) \times N_{cluster}^{-1}}
}
\item{imp}{Name or column number of the imputation variable.
}
  \item{group}{
%%     ~~Describe \code{sig.dif.bound} here~~
Optional: column number or name of one grouping variable. Note: in contrast to \code{repMean}, only one grouping variable can be specified.
}
  \item{trend}{
Optional: name or column number of the trend variable which contains the measurement time of the survey.
\code{repLmer} computes differences for all pairwise contrasts defined by trend variable levels. or three measurement
occasions, i.e. 2010, 2015, and 2020, contrasts (i.e. trends) are computed for 2010 vs. 2015, 2010 vs. 2020, and
2015 vs. 2020.
}
  \item{dependent}{
Name or column number of the dependent variable
}
  \item{formula.fixed}{
An R formula for fixed effects
}
  \item{formula.random}{
An R formula for random effects
}
  \item{doCheck}{
Logical: Check the data for consistency before analysis? If \code{TRUE} groups with insufficient data are excluded from
analysis to prevent subsequent functions from crashing.
}
  \item{na.rm}{
Logical: Should cases with missing values be dropped?
}
  \item{clusters}{
Variable name or column number of cluster variable.
}
  \item{verbose}{
Logical: Show analysis information on console?
}
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
\examples{
### load example data (long format)
data(lsa)
### use only the first nest, use only reading
btRead <- subset(lsa, nest==1 & domain=="reading")

\donttest{
### random intercept model with groups
mod1 <- repLmer(datL = btRead, ID = "idstud", wgt = "wgt", L1wgt="L1wgt", L2wgt="L2wgt",
        type = "jk2", PSU = "jkzone", repInd = "jkrep", imp = "imp",trend="year",
        group="country", dependent="score", formula.fixed = ~as.factor(sex)+mig,
        formula.random=~1, clusters="idclass")
res1 <- report(mod1)

### random slope without groups and without trend
mod2 <- repLmer(datL = subset(btRead, country=="countryA" & year== 2010),
        ID = "idstud", wgt = "wgt", L1wgt="L1wgt", L2wgt="L2wgt", type = "jk2",
        PSU = "jkzone", repInd = "jkrep", imp = "imp", dependent="score",
        formula.fixed = ~as.factor(sex)*mig, formula.random=~mig, clusters="idclass")
res2 <- report(mod2)
}}


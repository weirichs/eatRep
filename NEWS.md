# eatRep 0.15.1

* bugfix in the presentation of the regression results (reporting function for `repGlm()`)
* bugfix in the naming of the label column in the plain sheet of the reporting function `report2()`

# eatRep 0.15.0

* multicore computation is supported for all analyses using the survey package
* computation of group differences using 'group.differences.by' now also works for adjusted means
* for weighted analyses, the eatRep functions now remove zero weight cases prior to analyses to avoid counting them when determining sample size
* the eatRep functions now issue more informative alerts if missings occur on dependent/independent or jackknife variables. Missing occurrences concerned as crucial will be caught with an error message
* enhanced warning messages if number of nests/imputations differ between grouping variables
* enhanced warning messages if number of jkrep units differ between nests/imputations
* bugfix in the cross-level differences of group differences computation
* bugfix if levels of grouping variables contain leading and/or trailing spaces
* new function `report2` which summarizes the output of the four main functions a little tidier and provides an interface for eatPlot. The old reporting function `report()` is deprecated
* new function `pool.R2` computes pooled R^2 for multiple imputed and nested multiple imputed regression analyses according to Harel (2009)

# eatRep 0.14.7

* exemplary with variable labels stored as attributed
* add 'jkfac' argument to functions which call BIFIEsurvey functions
* trend analyses (with warning messages) now are possible when levels of grouping variables differ between assessment cycles
* bug fix in the internally used function which checks linking error object for consistency
* bug fix in repLmer when no weights are defined
* reporting function was slightly simplified (i.e., remove additional target formats)

# eatRep 0.14.6

* add function `repLmer()` for replication methods for linear multilevel models (wrapper to `BIFIEsurvey::BIFIE.twolevelreg()`)
* bug fix in internally used check function whether levels of grouping variable(s) are equal across assessment cycles (trend groups)
* bug fix in singularity treatment of `repGlm()`
* bug fix in `repQuantile()`
* add Ns (sample size) to the output of `repTable()` function
* add tests for `repQuantile()` and `repLmer()`

# eatRep 0.14.5

* compatibility with lavaan version 0.6-12

# eatRep 0.14.4

* bug fix in computation of adjusted weighted means
* enhance performance in recursive calls (e.g., trend estimation)

# eatRep 0.14.3

* add alternative output formats for the reporting function (argument `target`)
* bug fix in standard error computation of trend estimates
* bug fix in repQuantile when using 0% or 100% percentil
* enhance performance when using the BIFIEsurvey wrapper

# eatRep 0.14.0

* sorting of parameters and determination coefficient enhanced in reporting function for `repGlm()`
* add function `checkLEs` which checks consistency of linking errors with data which stem from an eatGADS data base
* trend estimation was broadened for more than two measurement occasions
* add class identifier in exemplary data

# eatRep 0.13.7

* add CR0 and CR2 methods for weighted effect coding with heterogeneous variances in `repMean()` and `repGlm()`, including a cluster argument
* add class ID variable to example data set
* bug fix: fixed inadequate column definition in the reporting function for `repQuantile()`
* adapt repQuantile function for rewritten function svyquantile from the survey package
* slightly revised vignette

# eatRep 0.13.6

* add asterisks for significance in the console output of the reporting function for `repGlm()`
* bug fix: fixed inadvertent warning message in function for group checks
* bug fix: call svyby with return.replicates = FALSE for quantile function

# eatRep 0.13.5

* added vignette
* added a balanced repeated replication example in the examples for `repMean()`
* added an internal check that persons are nested within groups

## Bug Fixes
* fixed `report` for cross differences when levels of grouping variable contained the variable name
* fixed missing entries in mode column in the output of  `repMean()`
* `eatRep::repMean(...)` now works as well as `library(eatRep); repMean(...)` 
* added tests for `repGlm()` function

# eatRep 0.13.4.9000

* Switch to Github Action for CI

# eatRep 0.13.4

* Initial CRAN release

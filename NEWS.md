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

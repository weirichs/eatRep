## Bug Fixes
* fixed inadvertent warning message in function for group checks

# eatRep 0.13.6

* add asterisks for significance in the console output of the reporting function for `repGlm()`

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

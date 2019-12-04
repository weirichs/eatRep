
# setwd("C:/Benjamin_Becker/02_Repositories/packages/eatRep")

context("cross level differences")

#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction.RData")
load("helper_SE_correction.RData")
#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction_complex.RData")
load("helper_SE_correction_complex.RData")
#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction_table.RData")
load("helper_SE_correction_table.RData")


test_that("Group variable creation for crossDiffs", {
  out <- suppressWarnings(report(means4T))
  cross_names <- unique(out[out$comparison == "crossDiff", "group"])
  crossgroup_names <- unique(out[out$comparison == "crossDiff_of_groupDiff", "group"])
  expect_true("female.vs.wholeGroup" %in% cross_names)
  expect_true("LandA_FALSE.vs.wholeGroup"  %in%  cross_names)
  expect_true("LandC_female_FALSE.vs.female_FALSE"  %in%  cross_names)
  expect_true("LandC_female_FALSE.vs.LandC_FALSE"  %in%  cross_names)  # fehlt bisher
  expect_true("LandC_female_FALSE.vs.LandC_female"  %in%  cross_names)
})


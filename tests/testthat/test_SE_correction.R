
# setwd("C:/Benjamin_Becker/02_Repositories/packages/eatRep")

context("SE correction for WEC and PISA")

#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction.RData")
load("helper_SE_correction.RData")
#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction_complex.RData")
load("helper_SE_correction_complex.RData")


test_that("Warnings for not supported crossDiffs", {
  warns <- capture_warnings(report(means3T))
  expect_equal(warns[[1]], "Standard error correction for 'crossDiff_of_groupDiff' is currently not supported.")
  expect_equal(warns[[2]], "Standard error correction for crossDifferences across multiple hierarchy levels is currently not supported.")
  #wec_out <- report(m_wec)
  #expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "se"], 1.954211, tolerance=1e-3)
  #expect_error(report(m_pisa), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
  # expect_error(report(m_wec), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
})

test_that("Temporary error", {
  #expect_error(report(m_wec), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
  #wec_out <- report(m_wec)
  #expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "se"], 1.954211, tolerance=1e-3)
  expect_error(report(m_pisa), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
  # expect_error(report(m_wec), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
})

test_that("wec 1 grouping variable, no trend", {
  wec_out <- report(m_wec)
  expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "se"], 1.954, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "male.vs.wholeGroup" & wec_out$parameter == "mean", "p"], 0, tolerance=1e-3)
})

test_that("wec 2 grouping variables, trend", {
  wec_out <- suppressWarnings(report(means3T))
  expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "se_2010"], 3.063, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "male.vs.wholeGroup" & wec_out$parameter == "mean", "se_2015"], 2.830, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "LandA.vs.wholeGroup" & wec_out$parameter == "mean", "p_2010"], 0.003, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "LandC.vs.wholeGroup" & wec_out$parameter == "mean", "p_2015"], 0.017, tolerance=1e-3)
})

test_that("wec 2 grouping variables, no trend, cross diff on higher level", {
  wec_out <- suppressWarnings(report(means3Tb))
  expect_equal(wec_out[wec_out$group == "LandA.vs.wholeGroup" & wec_out$parameter == "mean", "se"], 0.763, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "p"], 0.001, tolerance=1e-3)
  
  expect_equal(wec_out[wec_out$group == "LandA.vs.LandA_female" & wec_out$parameter == "mean", "se"], 3.344, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "LandA.vs.LandA_male" & wec_out$parameter == "mean", "p"], 0.002, tolerance=1e-3)
  
  expect_equal(wec_out[wec_out$group == "LandC.vs.LandC_female" & wec_out$parameter == "mean", "p"], 0.026, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "LandC.vs.LandC_male" & wec_out$parameter == "mean", "se"], 3.325, tolerance=1e-3)
  # expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "p_2010"], 0.003, tolerance=1e-3)
  # expect_equal(wec_out[wec_out$group == "LandC.vs.wholeGroup" & wec_out$parameter == "mean", "p_2015"], 0.017, tolerance=1e-3)
})

test_that("Old method still works", {
  old_rep <- report(m_old)
  expect_equal(nrow(old_rep), 10)
  expect_equal(length(which(old_rep$comparison == "crossDiff")), 4)
  expect_equal(old_rep[which(old_rep$comparison == "crossDiff"), "group"], c("female.vs.wholeGroup", "female.vs.wholeGroup", "male.vs.wholeGroup", "male.vs.wholeGroup"))
})


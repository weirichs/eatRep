
# setwd("C:/Benjamin_Becker/02_Repositories/packages/eatRep")

context("SE correction for WEC and PISA")

#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SEcorrection.RData")
load("helper_SE_correction.RData")

#str(m_old)
#str(m_pisa$pisa)
#str(m_wec$wec)


test_that("dummy", {
  expect_equal(1+1, 2)
})


test_that("Temporary error", {
  expect_error(report(m_wec), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
  expect_error(report(m_pisa), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
})


test_that("Old method still works", {
  old_rep <- report(m_old)
  expect_equal(nrow(old_rep), 10)
  expect_equal(length(which(old_rep$comparison == "crossDiff")), 4)
  expect_equal(old_rep[which(old_rep$comparison == "crossDiff"), "group"], c("female.vs.wholeGroup", "female.vs.wholeGroup", "male.vs.wholeGroup", "male.vs.wholeGroup"))
})

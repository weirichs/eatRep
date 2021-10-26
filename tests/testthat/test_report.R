

#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction.RData")
#load("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatRep/tests/testthat/helper_SE_correction.RData")
load("helper_SE_correction.RData")
#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction_wec_complex.RData")
#load("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatRep/tests/testthat/helper_SE_correction_wec_complex.RData")
load("helper_SE_correction_wec_complex.RData")

test_that("reporting", {
  out <- report(m_wec)
  expect_true(all(c("female", "male", "female.vs.wholeGroup", "male.vs.wholeGroup", "wholeGroup") %in% out$group))
  expect_true(all(c("crossDiff") %in% out$comparison))
  expect_true(all(c("mean", "sd") %in% out$parameter))
})

test_that("reporting with level variable name potential conflicts", {
  m_wec2 <- m_wec
  m_wec2$resT$noTrend$group <- gsub("female", "sexF", m_wec2$resT$noTrend$group)
  m_wec2$resT$noTrend$group <- gsub("male", "sexM", m_wec2$resT$noTrend$group)
  m_wec2$resT$noTrend$sex <- gsub("female", "sexF", m_wec2$resT$noTrend$sex)
  m_wec2$resT$noTrend$sex <- gsub("male", "sexM", m_wec2$resT$noTrend$sex)
  m_wec2$SE_correction[[1]]$resT$noTrend$parameter <- gsub("female", "sexF", m_wec2$SE_correction[[1]]$resT$noTrend$parameter)
  m_wec2$SE_correction[[1]]$resT$noTrend$parameter <- gsub("male", "sexM", m_wec2$SE_correction[[1]]$resT$noTrend$parameter)
  
  out <- report(m_wec2)
  out2 <- report(m_wec)
  expect_true(all(c("sexF", "sexM", "sexF.vs.wholeGroup", "sexM.vs.wholeGroup", "wholeGroup") %in% out$group))
  expect_true(all(c("crossDiff") %in% out$comparison))
  expect_true(all(c("mean", "sd") %in% out$parameter))
  expect_equal(out[-c(1, 6)], out2[-c(1, 6)])
})


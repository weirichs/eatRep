

#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction.RData")
#load("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatRep/tests/testthat/helper_SE_correction.RData")
load("helper_SE_correction.RData")
#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction_wec_complex.RData")
#load("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatRep/tests/testthat/helper_SE_correction_wec_complex.RData")
load("helper_SE_correction_wec_complex.RData")
#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction_pisa_complex.RData")
#load("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatRep/tests/testthat/helper_SE_correction_pisa_complex.RData")
load("helper_SE_correction_pisa_complex.RData")
#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction_table.RData")
#load("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatRep/tests/testthat/helper_SE_correction_table.RData")
load("helper_SE_correction_table.RData")
#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction_brr.RData")
#load("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatRep/tests/testthat/helper_SE_correction_brr.RData")
load("helper_SE_correction_brr.RData")
#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_SE_correction_others.RData")
#load("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatRep/tests/testthat/helper_SE_correction_others.RData")
load("helper_SE_correction_others.RData")
#load("c:/Benjamin_Becker/02_Repositories/packages/eatRep/tests/testthat/helper_different_modi.RData")
#load("c:/diskdrv/Winword/Psycho/IQB/Dropbox/R/eat/eatRep/tests/testthat/helper_different_modi.RData")
load("helper_different_modi.RData")


test_that("Warnings for not supported crossDiffs", {
  warns <- capture_warnings(report(means3T))
  expect_equal(warns[[1]], "Standard error correction for 'crossDiff_of_groupDiff' is currently not supported.")
  expect_equal(warns[[2]], "Standard error correction for crossDifferences across multiple hierarchy levels is currently not supported.")
  #wec_out <- report(m_wec)
  #expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "se"], 1.954211, tolerance=1e-3)
  #expect_error(report(m_pisa), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
  # expect_error(report(m_wec), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
})

#test_that("Temporary error", {
  #expect_error(report(m_wec), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
  #wec_out <- report(m_wec)
  #expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "se"], 1.954211, tolerance=1e-3)
  #expect_error(report(m_pisa), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
  # expect_error(report(m_wec), "SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
#})

test_that("wec 1 grouping variable, no trend", {
  wec_out <- report(m_wec)
  expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "se"], 1.975, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "male.vs.wholeGroup" & wec_out$parameter == "mean", "p"], 0, tolerance=1e-3)
})
test_that("pisa 1 grouping variable, no trend", {
  pisa_out <- report(m_pisa)
  expect_lt(pisa_out[pisa_out$group == "female.vs.wholeGroup" & pisa_out$parameter == "mean", "se"], 3.495)
  expect_lte(pisa_out[pisa_out$group == "male.vs.wholeGroup" & pisa_out$parameter == "mean", "p"], 0.05)
})

test_that("wec 2 grouping variables, trend", {
  wec_out <- suppressWarnings(report(means3T))
  expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "se_2010"], 3.063, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "male.vs.wholeGroup" & wec_out$parameter == "mean", "se_2015"], 2.850, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "countryA.vs.wholeGroup" & wec_out$parameter == "mean", "p_2010"], 0.002, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "countryC.vs.wholeGroup" & wec_out$parameter == "mean", "p_2015"], 0.017, tolerance=1e-3)
})

# means3T_old <- means3T_pisa
# means3T_old[["SE_correction"]] <- NULL
# report(means3T_old)[which(report(means3T_old)$comparison == "crossDiff"), ]
# report(means3T)[which(report(means3T)$comparison == "crossDiff"), ]
# report(means3T_pisa)[which(report(means3T_pisa)$comparison == "crossDiff"), ]

# means3T_pisa$resT$`2010`[which(means3T_pisa$resT$`2010`$comparison == "crossDiff"),]
test_that("pisa 2 grouping variables, trend", {
  pisa_out <- suppressWarnings(report(means3T_pisa))
  expect_equal(pisa_out[pisa_out$group == "female.vs.wholeGroup" & pisa_out$parameter == "mean", "se_2010"], 3.066)
  expect_equal(pisa_out[pisa_out$group == "male.vs.wholeGroup" & pisa_out$parameter == "mean", "se_2015"], 2.850)
  expect_equal(pisa_out[pisa_out$group == "countryA.vs.wholeGroup" & pisa_out$parameter == "mean", "p_2010"], 0.002)
  expect_equal(pisa_out[pisa_out$group == "countryC.vs.wholeGroup" & pisa_out$parameter == "mean", "p_2015"], 0.017)
})


test_that("wec 2 grouping variables, no trend, cross diff on higher level", {
  wec_out <- suppressWarnings(report(means3Tb))
  table(wec_out[which(wec_out$comparison == "crossDiff"), "group"])
  expect_equal(wec_out[wec_out$group == "countryA.vs.wholeGroup" & wec_out$parameter == "mean", "se"], 0.776, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "p"], 0.001, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "countryA_female.vs.countryA" & wec_out$parameter == "mean", "se"], 3.344, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "countryA_male.vs.countryA" & wec_out$parameter == "mean", "p"], 0.002, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "countryC_female.vs.countryC" & wec_out$parameter == "mean", "p"], 0.026, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "countryC_male.vs.countryC" & wec_out$parameter == "mean", "se"], 3.366, tolerance=1e-3)
})

test_that("pisa 2 grouping variables, no trend, cross diff on higher level", {
  #expect_error(suppressWarnings(report(means3Tb_pisa)), "PISA method for SE correction has not been fully implemented yet. Use crossDiffSE = 'old'.")
  pisa_out <- suppressWarnings(report(means3Tb_pisa))

  #means3Tb_pisa$SE_correction[[1]]$resT
  expect_equal(pisa_out[pisa_out$group == "countryA.vs.wholeGroup" & pisa_out$parameter == "mean", "se"], 0.776, tolerance=1e-3)
  expect_equal(pisa_out[pisa_out$group == "female.vs.wholeGroup" & pisa_out$parameter == "mean", "se"], 2.953, tolerance=1e-3)
  expect_equal(pisa_out[pisa_out$group == "countryA_female.vs.countryA" & pisa_out$parameter == "mean", "se"], 3.346, tolerance=1e-3)
  expect_equal(pisa_out[pisa_out$group == "countryA_male.vs.countryA" & pisa_out$parameter == "mean", "p"], 0.002, tolerance=1e-3)
  expect_equal(pisa_out[pisa_out$group == "countryC_female.vs.countryC" & pisa_out$parameter == "mean", "p"], 0.025, tolerance=1e-3)
  expect_equal(pisa_out[pisa_out$group == "countryC_male.vs.countryC" & pisa_out$parameter == "mean", "se"], 3.366, tolerance=1e-3)
})


test_that("wec 3 grouping variables, no trend, cross diff on higher level", {
  wec_out <- suppressWarnings(report(means4T))
  # lapply(means4T$SE_correction, function(x) x$refGrp)
  expect_equal(wec_out[wec_out$group == "countryA.vs.wholeGroup" & wec_out$parameter == "mean", "se"], 0.776, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "female.vs.wholeGroup" & wec_out$parameter == "mean", "p"], 0.001, tolerance=1e-3)
  
  #report(means4T$SE_correction[[4]])
  expect_equal(wec_out[wec_out$group == "countryA_female_TRUE.vs.countryA_female" & wec_out$parameter == "mean", "se"], 7.110, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "countryA_male_FALSE.vs.countryA_male" & wec_out$parameter == "mean", "p"], 0, tolerance=1e-3)
  
  #report(means4T$SE_correction[[18]])
  expect_equal(wec_out[wec_out$group == "countryA_female_TRUE.vs.female_TRUE" & wec_out$parameter == "mean", "se"], 1.454, tolerance=1e-3)
  expect_equal(wec_out[wec_out$group == "countryC_female_TRUE.vs.female_TRUE" & wec_out$parameter == "mean", "p"], 0.030, tolerance=1e-3)
})

test_that("Warning for different point estimates", {
  expect_warning(compare_point_estimates(old_est = 5, new_est = 2, param = "vgl"),
                 "Difference in point estimate of cross level difference for comparison vgl: 3")
  
  m_pisa2 <- m_pisa
  m_wec2 <- m_wec
  m_pisa2$SE_correction[[1]]$resT$noTrend[1, "value"] <- 0.5
  expect_warning(report(m_pisa2))
  
  m_wec2$SE_correction[[1]]$resT$noTrend[1, "value"] <- 0.5
  expect_warning(report(m_wec2))
})

test_that("Point estiamates for BRR are ok", {
  expect_silent(report(m_pisa_brr))
  expect_silent(report(m_wec_brr))
})

test_that("Old method still works", {
  old_rep <- report(m_old)
  expect_equal(nrow(old_rep), 10)
  expect_equal(length(which(old_rep$comparison == "crossDiff")), 4)
  expect_equal(old_rep[which(old_rep$comparison == "crossDiff"), "group"], c("female.vs.wholeGroup", "female.vs.wholeGroup", "male.vs.wholeGroup", "male.vs.wholeGroup"))
})

test_that("Brackets are supported in factor levels", {
  expect_silent(old_rep <- suppressWarnings(report(m_brack)))
})

test_that("SE for SD still unaffected", {
  old <- report(m_old)
  wec <- report(m_wec)
  expect_equal(old[old$parameter == "sd", ], wec[wec$parameter == "sd", ])
  
  old2 <- suppressWarnings(report(means4T))
  wec2 <- report(means4T_old)
  expect_equal(old2[old2$parameter == "sd", ], wec2[wec2$parameter == "sd", ])
})

test_that("Table and other functions not affected (protection against implementation from Sebastian)", {
  # internal check
  jk2.out <- m_table
  expect_equal(!is.null(jk2.out[["SE_correction"]]) && !is.null(jk2.out[["SE_correction"]][[1]]), FALSE)
  
  # reporting unaffected
  old_table <- report(m_table)
  expect_equal(nrow(old_table), 25)
  expect_equal(length(which(old_table$comparison == "crossDiff")), 10)
  expect_equal(old_table[which(old_table$comparison == "crossDiff"), "group"], c(rep("female.vs.wholeGroup", 5), rep("male.vs.wholeGroup", 5)))
})

test_that("Modus changed for WEC", {
  out <- report(different_modi)
  expect_equal(unique(out[which(out$comparison == "crossDiff" & out$parameter == "mean"), "modus"]), c("JK1.glm"))
  expect_equal(unique(out[which(out$comparison != "crossDiff" | out$parameter != "mean"), "modus"]), c("JK1.mean__BIFIEsurvey"))
})

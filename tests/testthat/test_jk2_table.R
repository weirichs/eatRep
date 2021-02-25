

rd     <- lsa[which(lsa[,"domain"] == "reading"),]
rdN1   <- rd[which(rd[,"nest"] == 1),]
rdN1y10<- rdN1[which(rdN1[,"year"] == 2010),]

test_that("No se_correction for jk2.table", {
  suppressMessages(tab1_old <- repTable(datL = rdN1y10, ID = "idstud", wgt = "wgt", imp="imp",
                       type = "JK2", PSU = "jkzone", repInd = "jkrep", groups = "country", group.splits = 0:1,
                       group.differences.by = "country", dependent = "comp", chiSquare = FALSE, cross.differences = TRUE,
                       verbose=FALSE, progress = FALSE, crossDiffSE = "old"))
  mess <- capture_messages(tab1_wec <- repTable(datL = rdN1y10, ID = "idstud", wgt = "wgt", imp="imp",
                       type = "JK2", PSU = "jkzone", repInd = "jkrep", groups = "country", group.splits = 0:1,
                       group.differences.by = "country", dependent = "comp", chiSquare = FALSE, cross.differences = TRUE,
                       crossDiffSE = "wec", verbose=FALSE, progress = FALSE))
  expect_equal(mess, "To date, only method 'old' is applicable for cross level differences in frequency tables.\n")
  expect_equal(tab1_old$SE_correction[[1]], NULL)
  expect_equal(tab1_wec$SE_correction[[1]], NULL)
})


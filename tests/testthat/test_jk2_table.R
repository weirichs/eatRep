

context("jk2 table")

#library(eatRep)
data("lsa")
rd     <- lsa[which(lsa[,"domain"] == "reading"),]
rdN1   <- rd[which(rd[,"nest"] == 1),]
rdN1y10<- rdN1[which(rdN1[,"year"] == 2010),]

txt <- capture.output ( tab1_old <- jk2.table(datL = rdN1y10, ID = "idstud", wgt = "wgt", imp="imp",
                       type = "JK2", PSU = "jkzone", repInd = "jkrep", groups = "country", group.splits = 0:1,
                       group.differences.by = "country", dependent = "comp", chiSquare = FALSE, cross.differences = TRUE))

txt <- capture.output ( tab1_wec <- jk2.table(datL = rdN1y10, ID = "idstud", wgt = "wgt", imp="imp",
                                              type = "JK2", PSU = "jkzone", repInd = "jkrep", groups = "country", group.splits = 0:1,
                                              group.differences.by = "country", dependent = "comp", chiSquare = FALSE, cross.differences = TRUE,
                                              crossDiffSE = "wec"))

test_that("No se_correction for jk2.table", {
  expect_equal(tab1_old$SE_correction[[1]], NULL)
  expect_equal(tab1_wec$SE_correction[[1]], NULL)
})


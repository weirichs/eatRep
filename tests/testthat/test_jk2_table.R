

context("jk2 table")

#library(eatRep)
data("lsa")
rd     <- lsa[which(lsa[,"domain"] == "reading"),]
rdN1   <- rd[which(rd[,"nest"] == 1),]
rdN1y10<- rdN1[which(rdN1[,"year"] == 2010),]

txt <- capture.output ( freq.tab1 <- jk2.table(datL = rdN1y10, ID = "idstud", wgt = "wgt", imp="imp",
                       type = "JK2", PSU = "jkzone", repInd = "jkrep", groups = "country", group.splits = 0:1,
                       group.differences.by = "country", dependent = "comp", chiSquare = FALSE, cross.differences = TRUE))

test_that("No se_correction for jk2.table", {
  expect_equal(freq.tab1$SE_correction, NULL)
})


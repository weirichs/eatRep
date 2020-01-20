
context("jk2 mean")

data("lsa")
rd     <- lsa[which(lsa[,"domain"] == "reading"),]
rd15 <- rd[rd$year == 2015, ]

txt <- capture.output ( m_withoutCross <- jk2.mean(datL = rd15, ID="idstud", type = "JK2", PSU = "jkzone", repInd = "jkrep",
              imp="imp", nest="nest", groups = c("sex"), group.splits = 0:1,
              cross.differences = FALSE, dependent = "score", na.rm=FALSE, doCheck=TRUE, linkErr = "leScore", crossDiffSE="old",
              engine = "BIFIEsurvey"))

txt2 <- capture.output ( m_oldCross <- jk2.mean(datL = rd15, ID="idstud", type = "JK2", PSU = "jkzone", repInd = "jkrep",
              imp="imp", nest="nest", groups = c("sex"), group.splits = 0:1,
              cross.differences = TRUE, dependent = "score", na.rm=FALSE, doCheck=TRUE, linkErr = "leScore", crossDiffSE="old",
              engine = "BIFIEsurvey"))

test_that("No cross differences", {
  expect_equal(m_withoutCross[["SE_correction"]], NULL)
  expect_false("SE_correction" %in% names(m_withoutCross))
})


test_that("Old cross differences", {
  expect_equal(class(m_oldCross[["SE_correction"]]), c("old", "list"))
  expect_equal(m_oldCross[["SE_correction"]][[1]], NULL)
})

rd15$sex_logic <- as.logical(as.numeric(rd15$sex) - 1)

test_that("error for two logical grouping variables", {
  expect_error(capture.output(jk2.mean(datL = rd15, ID="idstud", type = "JK2", PSU = "jkzone", repInd = "jkrep",
           imp="imp", nest="nest", groups = c("sex_logic", "mig"), group.splits = 0:1,
           cross.differences = FALSE, dependent = "score", na.rm=FALSE, doCheck=TRUE, linkErr = "leScore", crossDiffSE="old",
           engine = "BIFIEsurvey")), "Factor levels of grouping variables are not disjunct.")
})


### PISA method
test_that("PISA runs through", {
  expect_silent(suppressWarnings(suppressMessages(txt2 <- capture.output ( m_oldCross <- jk2.mean(datL = rd15, ID="idstud", type = "JK2", PSU = "jkzone", repInd = "jkrep",
                                                                imp="imp", nest="nest", groups = c("sex"), group.splits = 0:1,
                                                                cross.differences = TRUE, dependent = "score", na.rm=FALSE, doCheck=TRUE, linkErr = "leScore", crossDiffSE="rep")))))
})


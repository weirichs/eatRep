bt         <- lsa[which(lsa[,"nest"] == 1),]
bt2010     <- bt[which(bt[,"year"] == 2010),]
bt2010read <- bt2010[which(bt2010[,"domain"] == "reading"),]


test_that("repGlm", {
  suppressWarnings(txt <- capture_output(mod1 <- repGlm(datL = bt2010read, ID = "idstud", wgt = "wgt", type = "jk2",
                 PSU = "jkzone", repInd = "jkrep", imp = "imp", groups = "country",
                 formula = score~sex, family ="gaussian")))
  txt2 <- capture_output(res1 <- report(mod1, printGlm = TRUE))
  expect_equal(res1[which(res1[,"parameter"] == "(Intercept)"),"est"], c(508.406, 502.607, 526.287))
  expect_equal(dim(res1), c(18, 9))
  expect_equal(unique(res1$group), c("countryA", "countryB", "countryC"))
  expect_equal(unique(res1$parameter), c("(Intercept)", "Ncases", "sexmale", "Nvalid" , "R2", "R2nagel"  ))
})




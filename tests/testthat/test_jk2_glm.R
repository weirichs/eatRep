bt         <- lsa[which(lsa[,"nest"] == 1),]
bt2010     <- bt[which(bt[,"year"] == 2010),]
bt2010read <- bt2010[which(bt2010[,"domain"] == "reading"),]


test_that("repGlm", {
  suppressWarnings(txt <- capture_output(mod1 <- repGlm(datL = bt2010read, ID = "idstud", wgt = "wgt", type = "jk2",
                 PSU = "jkzone", repInd = "jkrep", imp = "imp", groups = "country",
                 formula = score~sex, family ="gaussian")))
  txt2 <- capture_output(res1 <- report2(mod1, printGlm = TRUE)[["plain"]])
  expect_equal(res1[which(res1[,"parameter"] == "(Intercept)"),"est"], c(508.406, 502.607, 526.287))
  expect_equal(unique(res1$parameter), c("(Intercept)", "Ncases", "sexmale", "Nvalid" , "R2"  ))
})


### Example 4: weighted effect coding to estimate whether a specific country's mean
### differs from the overall mean (whereas the overall population is a composite of
### all countries). The procedure adapts the weighted effect coding procedures
### described in te Grotenhuis (2017) for multiple imputation and replication methods.
mod4 <- repGlm(datL = bt2010read, ID = "idstud", wgt = "wgt", type = "jk2", PSU = "jkzone", repInd = "jkrep", imp = "imp", formula = score~country, useWec=TRUE)
suppressWarnings(res4 <- report(mod4, printGlm = FALSE))
res4A<- report2(mod4, printGlm = FALSE)


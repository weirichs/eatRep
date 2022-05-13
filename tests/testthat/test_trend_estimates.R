# test of correct trend estimates and trend standard errors
# die Vorlagenwerte kommen aus eatRep 0.8.10
data(lsa)
rd     <- lsa[intersect(intersect(intersect(which(lsa[,"domain"] == "reading"),which(lsa[,"year"] %in% c(2010, 2015))), which(lsa[,"nest"] == 1)),which(lsa[,"country"] == "countryA")),]
res    <- repMean(datL = rd, ID="idstud", wgt="wgt", type = "JK2", PSU = "jkzone", repInd = "jkrep",imp="imp", dependent = "score", na.rm=FALSE, doCheck=TRUE, trend = "year", linkErr = "leScore")
res1   <- report(res)

test_that("reporting", {
    expect_equal(res1[which(res1[,"parameter"] == "mean"),"est_trend_2010.vs.2015"], -3.590147, tolerance=1e-3)
    expect_equal(res1[which(res1[,"parameter"] == "mean"),"se_2010"], 3.514924, tolerance=1e-3)
    expect_equal(res1[which(res1[,"parameter"] == "mean"),"se_2015"], 3.816191, tolerance=1e-3)
    expect_equal(res1[which(res1[,"parameter"] == "mean"),"se_trend_2010.vs.2015"], 5.253851, tolerance=1e-3)
})



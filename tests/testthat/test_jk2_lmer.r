# replicate BIFIEsurvey example analyses
data(data.bifie01, package="BIFIEsurvey")

# create dataset with replicate weights and plausible values
bdat1 <- BIFIEsurvey::BIFIE.data.jack( data=data.bifie01, jktype="JK_TIMSS", jkzone="JKCZONE",
            jkrep="JKCREP", wgt="TOTWGT", pv_vars=c("ASMMAT","ASSSCI") )


test_that("repLmer checks", {
    skip_on_cran()
    mod2c <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bdat1, dep="ASMMAT",  formula.fixed=~  female +  ASBG06A, formula.random=~ ASBG06A,
                idcluster="idschool", wgtlevel2="SCHWGT", maxiter=500, se=TRUE)
    datL  <- na.omit(eatTools::wideToLong(datWide = data.bifie01,noImp = c("idschool", "idstud", "JKCZONE", "JKCREP", "TOTWGT", "SCHWGT", "female", "ASBG06A"),
                imp = list(mat = paste0("ASMMAT0",1:5), sci = paste0("ASSSCI0", 1:5)) ))
    mod2d <- repLmer(datL=datL, ID="idstud", wgt = "TOTWGT", L2wgt="SCHWGT", type = "JK2",  PSU = "JKCZONE", repInd = "JKCREP", imp="imp",
                dependent="mat", formula.fixed=~female + ASBG06A, formula.random=~ASBG06A, doCheck = TRUE, na.rm = FALSE, clusters="idschool", verbose = TRUE)
    res2d <- report(mod2d)
    expect_equal(round(subset(res2d, parameter=="(Intercept)")[,"est"],digits=2), round(subset(mod2c[["stat"]], parameter=="beta_(Intercept)")[,"est"], digits=2))
    expect_equal(round(subset(res2d, parameter=="ASBG06A")[,"est"],digits=1), round(subset(mod2c[["stat"]], parameter=="beta_ASBG06A")[,"est"], digits=1))
    expect_equal(round(subset(res2d, parameter=="female")[,"est"],digits=1), round(subset(mod2c[["stat"]], parameter=="beta_female")[,"est"], digits=1))
    expect_equal(round(subset(res2d, parameter=="ICC_Cond")[,"est"],digits=3), round(subset(mod2c[["stat"]], parameter=="ICC_Cond")[,"est"], digits=3))
    expect_equal(round(subset(res2d, parameter=="ICC_Uncond")[,"est"],digits=3), round(subset(mod2c[["stat"]], parameter=="ICC_Uncond")[,"est"], digits=3))
    expect_equal(round(subset(res2d, parameter=="R2_Lev2")[,"est"],digits=3), round(subset(mod2c[["stat"]], parameter=="R2_Lev2")[,"est"], digits=3))
    expect_equal(round(subset(res2d, parameter=="R2_Lev1")[,"est"],digits=3), round(subset(mod2c[["stat"]], parameter=="R2_Lev1")[,"est"], digits=3))
})
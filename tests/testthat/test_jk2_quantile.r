rd     <- lsa[which(lsa[,"domain"] == "reading"),]

test_that("complete quantile checks", {
    skip_on_cran()
    ### Varianten durchpermutieren, dauert ca. 30 Sek.
    typ    <- do.call("rbind", lapply(c("none", "jk2"), FUN = function (typ) {
              if(typ == "none") {jkz <- NULL; jkr <- NULL } else {jkz <- "jkzone"; jkr <- "jkrep"}
              imp <- do.call("rbind", lapply(c("noImp", "imp", "nest"), FUN = function (i) {
                     if ( i == "nest") {dat <- rd; impV <- "imp"; nestV <- "nest"}
                     if ( i == "imp")   {dat <- subset(rd, nest == 1); impV <- "imp"; nestV <- NULL}
                     if ( i == "noImp")  {dat <- subset(rd, imp==1 & nest == 1); impV <- NULL; nestV <- NULL}
                     wgt <- do.call(plyr::rbind.fill, lapply(c("noWgt", "wgt"), FUN = function (w ) {
                            if ( w == "noWgt") {w1 <- NULL} else { w1 <- "wgt"}
                            perzent   <- repQuantile(datL = dat, ID = "idstud", wgt = w1, type = typ,  PSU = jkz, repInd = jkr,
                                         imp = impV, nest=nestV,  groups = "country", group.splits = c(0:1), dependent = "score",
                                         probs = c(0, 0.05,0.9,1) , trend = "year")
                            res       <- report2(perzent, add = list(domain = "reading", weights = w, imp=i))[["plain"]]
                            return(res)}))                                      ### perzentile 0 und 1 sind kritische Grenzen
                     return(wgt)}))
              return(imp)}))

    # dimension von rueckgabeobjekt
    expect_true(inherits(typ[,"parameter"], "character"))

    # parameterspalte numerisch machen
    typ[,"parameter"] <- as.numeric(typ[,"parameter"])

    # 2010
    typ2010 <- reshape2::dcast(subset(typ,comparison=="none" & year=="2010"), country+weights+imp+parameter~modus, value.var="est")
    sub2010 <- subset(typ2010, parameter>0 & parameter < 1)
    expect_true(max(abs(sub2010[,"CONV.quantile"] - sub2010[,"JK2.quantile__survey"])) < 0.75)
    expect_true(max(abs(typ2010[,"CONV.quantile"] - typ2010[,"JK2.quantile__survey"])) < 6)

    # 2015
    typ2015 <- reshape2::dcast(subset(typ,comparison=="none" & year=="2015"), country+weights+imp+parameter~modus, value.var="est")
    sub2015 <- subset(typ2015, parameter>0 & parameter < 1)
    expect_true(max(abs(sub2015[,"CONV.quantile"] - sub2015[,"JK2.quantile__survey"])) < 1)
    expect_true(max(abs(typ2015[,"CONV.quantile"] - typ2015[,"JK2.quantile__survey"])) < 5)

    # trend
    typTrend <- reshape2::dcast(subset(typ,comparison=="trend"), country+weights+imp+parameter~modus, value.var="est")
    subTrend <- subset(typTrend, parameter>0 & parameter < 1)
    expect_true(max(abs(subTrend[,"CONV.quantile"] - subTrend[,"JK2.quantile__survey"])) < 1.25)
    expect_true(max(abs(typTrend[,"CONV.quantile"] - typTrend[,"JK2.quantile__survey"])) < 2.8)
})
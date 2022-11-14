rd     <- lsa[which(lsa[,"domain"] == "reading"),]
rdN1   <- rd[which(rd[,"nest"] == 1),]
rdN1y10<- rdN1[which(rdN1[,"year"] == 2010),]

test_that("complete WEC checks", {
    skip_on_cran()
    ### alle Varianten durchpermutieren, insgesamt 8
    mods   <- do.call("rbind", lapply(c("lm", "lavaan"), FUN = function (eng) {
              m1 <- do.call("rbind", lapply(c(TRUE, FALSE), FUN = function ( het) {
                    m2 <- do.call("rbind", lapply(c("rep", "norep"), FUN = function (re) {
                          if ( re == "rep") {
                              type <- "JK2"; PSU <- "jkzone"; repInd <- "jkrep"
                          }  else  {
                              type <- "none"; PSU <- NULL; repInd <- NULL
                          }
                          means1a<- repMean(datL = rdN1y10, ID="idstud", wgt="wgt", type = type, PSU=PSU, repInd=repInd, imp="imp", groups = "country", group.splits = 0:1, group.differences.by = "country", cross.differences = TRUE, dependent = "score", na.rm=FALSE, doCheck=TRUE, hetero=het, crossDiffSE.engine= eng)
                          res    <- report(means1a, exclude="var", add = list(crossDiffSE.engine = eng, hetero = as.character(het)))
                          return(res)}))
                    return(m2)}))
              return(m1)}))
    expect_equal(dim(mods), c(200,12))
    wide1  <- reshape2::dcast(mods[grep("mean", mods[,"modus"]),], group+hetero+parameter~modus+crossDiffSE.engine, value.var="est")
    wide1[,"abw"] <- unlist(plyr::alply(wide1, .margins = 1, .fun = function (z) {sd(z[,c("CONV.mean_lavaan", "CONV.mean_lm", "JK2.mean__survey_lavaan", "JK2.mean__survey_lm")])}))
    expect_true(all(wide1[which(wide1[,"parameter"] != "sd"),"abw"] == 0))
    abw    <- wide1[which(wide1[,"parameter"] == "sd"),"abw"]
    expect_true(all(abs(abw) < 0.025))
})
    
          
          





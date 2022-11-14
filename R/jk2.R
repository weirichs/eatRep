createCall <- function ( hetero, allNam, formula) {
         part1 <- ifelse(hetero, yes = "estimatr::lm_robust(", no = "lm(")
         part2 <- "formula, data = dat.i"
         if ( is.null(allNam[["wgt"]]) && is.null(allNam[["clusters"]]) && !hetero) {
              part3 <- paste0(part1, part2, ")")
              return(part3)
         }
         if ( !is.null(allNam[["wgt"]])) {
                   part2 <- paste0(part2, ", weights = ",allNam[["wgt"]] )
         }
         if ( hetero) {
              part2 <- paste0(part2, ", se_type = se_type")
         }
         if(!is.null(allNam[["clusters"]])) {
              part2 <- paste0(part2, ", clusters = ",allNam[["clusters"]])
         }
         part3 <- paste0(part1, part2, ")")
         return(part3)}

identify_UV_AV <- function ( toCall, dependent, formula, formula.random, formula.fixed, glmerFormula)  {
          independent <- NULL                                                   
          if(toCall == "glm") {                                                 
             dependent  <- as.character(formula)[2]
             independent<- all.vars(formula[-2])
          }
          if(toCall == "lmer") {
             independent<- unique(c(all.vars(formula.random), all.vars(formula.fixed)))
          }
          if(toCall == "glmer") {
             dependent  <- as.character(glmerFormula)[2]
             independent<- all.vars(glmerFormula[-2])
          }
          return(list(uv=independent, av=dependent))}

extractFactorVarsFromFormula <- function ( formula) {
         r <- Reduce(paste0, deparse(formula))
         r2<- methods::el(regmatches(r, gregexpr("(?<=factor\\().*?(?=\\))", r, perl=TRUE)))
         return(r2)}

chooseSeType <- function ( se_type, clusters) {
       if ( length(se_type) > 1) {                                              
            if ( is.null(clusters)) { se_type <- "HC3" } else { se_type <- "CR2"}
       }
       se_type <- match.arg(arg = se_type, choices = c("HC3", "HC0", "HC1", "HC2", "CR0", "CR2"))
       return(se_type)}

generateRandomJk1Zones <- function (datL, unit, nZones, name = "randomCluster") {
       datL  <- eatTools::makeDataFrame(datL, name = "data")
       stopifnot(length(unit)==1)
       allVar<- list(ID = unit)
       allNam<- eatTools::existsBackgroundVariables(dat = datL, variable=unlist(allVar), warnIfMissing = FALSE)
       if ( "randomCluster" %in% colnames(datL)) {stop("Name '",name,"' already exists in data. Please choose an alternative name.")}
       if ( nZones >= length(unique(datL[,allNam])) ) { stop("Number of zones must not exceed number of units.")}
       if ( nZones >= length(unique(datL[,allNam])) / 5 ) {warning("Number of zones (",nZones,") is large compared to the number of distinct units (",length(unique(datL[,allNam])),").")}
       reps  <- length(unique(datL[,allNam])) / nZones
       zones <- rep(1:nZones, times = ceiling(reps))
       zones <- data.frame ( ID = sample(unique(datL[,allNam]), size = length(unique(datL[,allNam])), replace=FALSE), zone = zones[1:length(unique(datL[,allNam]))], stringsAsFactors = FALSE)
       colnames(zones)[2] <- name
       mdat  <- merge(datL, zones, by.x = allNam, by.y = "ID", all = TRUE)
       return(mdat)}

repMean <- function(datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"), PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL,
            group.splits = length(groups), group.differences.by = NULL, cross.differences = FALSE, crossDiffSE = c("wec", "rep","old"), adjust = NULL, useEffectLiteR = FALSE, nBoot = 100,
            group.delimiter = "_", trend = NULL, linkErr = NULL, dependent, na.rm = FALSE, doCheck = TRUE, engine = c("survey", "BIFIEsurvey"), scale = 1, rscales = 1, mse=TRUE, rho=NULL, hetero=TRUE, se_type = c("HC3", "HC0", "HC1", "HC2", "CR0", "CR2"),
            clusters =NULL, crossDiffSE.engine= c("lavaan", "lm"), stochasticGroupSizes = FALSE, verbose = TRUE, progress = TRUE) {
            depO <- attr(datL, "depOri")
            fc   <- attr(datL, "fc")
            crossDiffSE.engine <- match.arg(crossDiffSE.engine)
            cdse<- match.arg(arg = crossDiffSE, choices = c("wec", "rep","old"))
            if (!is.null(adjust)) {
                if(is.list(cross.differences) || isTRUE(cross.differences)) {
                     if ( cdse != "old") {
                         warning("To date, for adjusted means, cross-level differences can only be computed with method 'old'. Set 'crossDiffSE' to 'old'.")
                         cdse <- "old"
                     }
                }
            }
            type    <- car::recode(match.arg(arg = toupper(type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'")
            se_type <- chooseSeType(se_type, clusters)
            datL    <- eatTools::makeDataFrame ( datL, name = "datL")
            if ( is.null ( attr(datL, "modus"))) {
                  modus <- identifyMode ( name = "mean", type = car::recode(match.arg(arg = toupper(type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'"))
            }  else  {
                  modus <- attr(datL, "modus")
            }
            ret <- eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "mean",
                   nest = nest, imp = imp, groups = groups, group.splits = group.splits, group.differences.by = group.differences.by,
                   cross.differences = cross.differences, adjust=adjust, useEffectLiteR = useEffectLiteR, trend = trend, linkErr = linkErr, dependent = dependent, group.delimiter=group.delimiter, na.rm=na.rm, doCheck=doCheck, modus = modus, engine=engine,
                   scale = scale, rscales = rscales, mse=mse, rho=rho, crossDiffSE.engine= crossDiffSE.engine, stochasticGroupSizes=stochasticGroupSizes, verbose=verbose, progress=progress, clusters=clusters, fc=fc, depOri=depO)
            if ( isFALSE(ret[["allNam"]][["cross.differences"]])) {return(ret)}
            if ( (is.list(cross.differences) || cross.differences == TRUE) && cdse == "old" ) {
                 ret[["SE_correction"]] <- list(NULL)
                 class(ret[["SE_correction"]]) <- c("old", "list")
                 return(ret)
            }
            if ( is.list(cross.differences) || isTRUE(cross.differences) ) {
                  toAppl<- superSplitter(group = ret[["allNam"]][["group"]], group.splits = group.splits, group.differences.by = ret[["allNam"]][["group.differences.by"]], group.delimiter = group.delimiter , dependent=ret[["allNam"]][["dependent"]] )
                  stopifnot(length ( toAppl ) > 1)
                  vgl <- do.call("rbind", lapply(1:length(toAppl), FUN = function ( y ) {
                         gdb <- attr(toAppl[[y]], "group.differences.by")
                         if ( is.null(gdb)) {gdb <- NA}
                         res <-  data.frame ( analysis.number = y, hierarchy.level = length(toAppl[[y]]), groups.divided.by = paste(toAppl[[y]], collapse=" + "), group.differences.by = gdb)
                         return(res)}))
                  ana <- lapply ( combinat::combn(x=vgl[,"analysis.number"], m=2, simplify=FALSE), FUN = function ( a ) {
                         t1 <- abs(diff(vgl[ match(a, vgl[,"analysis.number"]),"hierarchy.level"])) == 1
                         t2 <- TRUE                                             
                         if ( vgl[min(a),"hierarchy.level"] != 0 ) {            
                              nam<- lapply(sort(a), FUN = function ( z ) {      
                                    n1 <- unlist(strsplit(as.character(vgl[z,"groups.divided.by"]), split=" |\\+"))
                                    n1 <- n1[which(nchar(n1)>0)]
                                    return(n1)})
                              if ( length(setdiff(nam[[2]], nam[[1]])) != 1) { t2 <- FALSE }
                         }
                         t3 <- TRUE
                         if (  is.list(cross.differences) ) {
                             if ( sum(unlist(lapply(cross.differences, FUN = function (v1) { all(sort(vgl[ match(a, vgl[,"analysis.number"]),"hierarchy.level"]) == sort(v1))}))) == 0) {t3 <- FALSE}
                         }
                         if(isTRUE(t1) && isTRUE(t2) && isTRUE(t3)) {return(a)} else {return(NULL)} })
                  ana <- ana[which(unlist(lapply(ana, FUN = function (l) {!is.null(l)}))==TRUE)]
                  if ( sum(abs(unlist(lapply(cross.differences, diff))) > 1) > 0 ) {
                      warning("Computation of cross level differences using '",cdse,"' method is only possible for differences according to adjacent levels. Non-adjacent levels will be ignored.")
                  }  else {
                      if ( verbose ) {
                            if (cdse != "wec" || is.null(ret[["allNam"]][["PSU"]])) {
                                cat(paste0("Compute cross level differences using '",cdse,"' method.\n"))
                            }  else  {
                                cat(paste0("Compute cross level differences using '",cdse,"' method. Assume ",car::recode(hetero, "TRUE='heteroscedastic'; FALSE='homoscedastic'")," variances.\n"))
                            }
                      }
                  }
                  spl <- lapply(ana, FUN = function ( a ) {
                         vgl <- vgl[a,]
                         if (vgl[1,"hierarchy.level"] != 0) {                   
                             nam <- unlist(strsplit(as.character(vgl[1,"groups.divided.by"]), split=" |\\+"))
                             nam <- nam[which(nchar(nam)>0)]
                             dat <- by(datL, INDICES = datL[,nam], FUN = function ( d ) {return(d)})
                             grp <- setdiff(unlist(strsplit(as.character(vgl[2,"groups.divided.by"]), split=" |\\+")) , nam)
                             grp <- grp[which(nchar(grp)>0)]
                         }  else  {                                             
                             dat <- list(datL)
                             grp <- as.character(vgl[2,"groups.divided.by"])
                         }
                         stopifnot(length(grp)==1)
                         cld <- lapply(dat, FUN = function ( d ) {              
                                if ( is.null(ret[["allNam"]][["wgt"]]) || !ret[["allNam"]][["wgt"]] %in% colnames(d)) {
                                     stopifnot(is.null(wgt))                    
                                     message("   '",cdse,"' method: Assume equally weighted cases.")
                                     gew <- NULL                                
                                } else {                                        
                                     gew <- ret[["allNam"]][["wgt"]]
                                }                                               
                                if (is.null(nest)) {ne <- NULL} else {ne <- ret[["allNam"]][["nest"]]}
                                if (is.null(imp))  {im <- NULL} else {im <- ret[["allNam"]][["imp"]]}
                                if ( vgl[1,"hierarchy.level"] != 0) {           
                                     rg <- eatTools::facToChar(d[1,nam, drop=FALSE])
                                     rg <- do.call("rbind", lapply(names(rg), FUN = function (r){data.frame ( groupName = r, groupValue = gsub("\\.", "", gsub("_", "",as.character(rg[1,r]))), stringsAsFactors = FALSE) }))
                                }  else  {
                                     rg <- NULL
                                }
                                if ( isTRUE(hetero)) {
                                     if ( cdse == "rep"){
                                         warning("Method 'rep' is not yet adapted for heterogeneous variances. Results might not be trustworthy.")
                                     }
                                }  else  {
                                     if(!is.null(ret[["allNam"]][["clusters"]])) {
                                         stop("If clusters are specified, 'hetero' must be TRUE.")
                                     }
                                }
                                if ( cdse == "wec" ) {                          
                                     if ( !inherits(d[,grp], "factor") ) {
                                         warning("Group variable '",grp,"' must be of class 'factor' for '",cdse,"'. Change class of '",grp,"' from '",class(d[,grp]),"' to 'factor'.")
                                         d[,grp] <- as.factor(d[,grp])
                                     }
                                     attr(d, "wrapperForWec") <- TRUE
                                     if ( isTRUE(stochasticGroupSizes) && ( !is.null(ret[["allNam"]][["PSU"]]) || !is.null(ret[["allNam"]][["repWgt"]]) ) ) {
                                          message("To date, stochastic group sizes cannot be used in combination with replication mehods. Switch to fixed group sizes.")
                                          stochasticGroupSizes <- FALSE
                                     }
                                     if ( isTRUE(stochasticGroupSizes) &&  crossDiffSE.engine == "lm" ) {
                                          message("To date, stochastic group sizes cannot be computed with 'lm' engine. Switch to crossDiffSE.engine == 'lavaan'.")
                                          crossDiffSE.engine  <- "lavaan"
                                     }
                                     b <- repGlm(datL=d, ID=ret[["allNam"]][["ID"]], wgt = gew, type = type, PSU = ret[["allNam"]][["PSU"]], repInd = ret[["allNam"]][["repInd"]],
                                                  repWgt = ret[["allNam"]][["repWgt"]], nest=ne, imp=im, trend = trend,
                                                  formula = as.formula(paste0(ret[["allNam"]][["dependent"]] , " ~ ", grp)), doCheck = doCheck, na.rm = na.rm, useWec = TRUE,
                                                  scale = scale, rscales = rscales, mse=mse, rho=rho, hetero=hetero, se_type=se_type , crossDiffSE.engine= crossDiffSE.engine, stochasticGroupSizes=stochasticGroupSizes,
                                                  verbose=verbose, progress=progress, clusters=clusters)
                                }  else  {                                      
                                     b <- eatRep(datL =d, ID=ret[["allNam"]][["ID"]], wgt = gew, type=type, PSU = ret[["allNam"]][["PSU"]], repInd = ret[["allNam"]][["repInd"]], toCall = "cov",nBoot=nBoot,
                                          nest=ne, imp=im, groups = grp, refGrp = rg, trend = trend, dependent = ret[["allNam"]][["dependent"]], na.rm=na.rm, doCheck=FALSE, engine="survey", modus=modus,scale = scale,
                                          rscales = rscales, mse=mse, rho=rho, reihenfolge = ret[["allNam"]][["group"]], verbose=verbose, progress=progress, clusters=clusters, fc=fc)
                                }
                                b[["vgl"]]      <- vgl
                                b[["focGrp"]]   <- grp                          
                                b[["refGrp"]]   <- "all"                        
                                if(!is.null(rg)) {b[["refGrp"]] <- rg}
                         return(b)})
                  return(cld)})
                  spl <- unlist(spl, recursive=FALSE)
                  class(spl) <- c( car::recode(cdse, "'wec'='wec_se_correction'; 'rep'='rep_se_correction'"), "list")
                  ret[["SE_correction"]] <- spl
            }
            return(ret)}

repTable<- function(datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"),
            PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL, group.splits = length(groups), group.differences.by = NULL, cross.differences = FALSE, crossDiffSE = c("wec", "rep","old"),
            nBoot = 100, chiSquare = FALSE, correct = TRUE, group.delimiter = "_", trend = NULL, linkErr = NULL, dependent , separate.missing.indicator = FALSE,na.rm=FALSE, expected.values = NULL, doCheck = TRUE, forceTable = FALSE,
            engine = c("survey", "BIFIEsurvey"), scale = 1, rscales = 1, mse=TRUE, rho=NULL, verbose = TRUE, progress = TRUE ) {
            crossDiffSE <- "old"                                                
            if(isFALSE(cross.differences) == FALSE) {message("To date, only method 'old' is applicable for cross level differences in frequency tables.")}
            modus <- identifyMode ( name = "table", type = car::recode(match.arg(arg = toupper(type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'"))
            datL  <- eatTools::makeDataFrame ( datL)
            chk1  <- eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "table",
                     nest = nest, imp = imp, groups = groups, group.splits = group.splits, group.differences.by = group.differences.by, cross.differences=cross.differences, correct = correct,
                     trend = trend, linkErr = linkErr, dependent = dependent, group.delimiter=group.delimiter, separate.missing.indicator=separate.missing.indicator,
                     expected.values=expected.values, na.rm=na.rm, doCheck=FALSE, onlyCheck= TRUE, modus = modus, engine=engine, scale = scale, rscales = rscales, mse=mse, rho=rho, verbose=verbose, progress=progress, clusters=NULL, fc="repTable")
            if ( length(unique(datL[,chk1[["dependent"]]])) == 2 && isTRUE(all(sort(unique(datL[,chk1[["dependent"]]])) == 0:1)) && isFALSE(forceTable)) {
                 attr(datL, "modus") <- modus
                 attr(datL,"fc") <- "repTable"
                 ret <- repMean ( datL = datL, ID=chk1[["ID"]], wgt=chk1[["wgt"]], type = type, PSU = chk1[["PSU"]], repInd = chk1[["repInd"]], repWgt = repWgt,
                        nest = chk1[["nest"]], imp = chk1[["imp"]], groups = groups, group.splits = group.splits, group.differences.by=group.differences.by, cross.differences=cross.differences,
                        crossDiffSE = crossDiffSE, nBoot = nBoot, group.delimiter =group.delimiter, trend = chk1[["trend"]], linkErr = chk1[["linkErr"]], dependent = chk1[["dependent"]],
                        na.rm=na.rm,doCheck = doCheck, engine=engine, scale = scale, rscales = rscales, mse=mse, rho=rho, verbose=verbose, progress=progress, clusters = NULL)
                 return(ret)
            }  else  {
                 if ( !is.null(group.differences.by) && isFALSE(chiSquare)) {
                    chk <- eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "table",
                           nest = nest, imp = imp, groups = groups, group.splits = group.splits, group.differences.by = group.differences.by, cross.differences=cross.differences, correct = correct,
                           trend = trend, linkErr = linkErr, dependent = dependent, group.delimiter=group.delimiter, separate.missing.indicator=separate.missing.indicator,
                           expected.values=expected.values, na.rm=na.rm, doCheck=doCheck, onlyCheck= TRUE, modus = modus, engine=engine, scale = scale, rscales = rscales,
                           mse=mse, rho=rho, verbose=verbose, progress=progress, clusters = NULL, fc="repTable")
                    isNa<- which ( is.na ( datL[, chk[["dependent"]] ] ))
                    if ( length ( isNa ) > 0 ) {
                         warning("Warning: Found ",length(isNa)," missing values in dependent variable '",chk[["dependent"]],"'.")
                         if ( isTRUE(separate.missing.indicator) ) {
                              stopifnot ( length( intersect ( "missing" , names(table(datL[, chk[["dependent"]] ])) )) == 0 )
                              if(inherits(datL[, chk[["dependent"]] ], "factor")){
                                  levOld <- levels(datL[, chk[["dependent"]] ]) ### dat[which(is.na(dat[,"var"])) ,"var"] <- "missing" nicht
                                  datL[, chk[["dependent"]] ] <- as.character(datL[, chk[["dependent"]] ])
                                  datL[isNa, chk[["dependent"]] ] <- "missing"
                                  datL[, chk[["dependent"]] ] <- factor(datL[, chk[["dependent"]] ], levels=c(levOld, "missing"))
                              }  else  {
                                  datL[isNa, chk[["dependent"]] ] <- "missing"
                              }
                         }  else  {
                              if ( isFALSE(na.rm ) ) { stop("If no separate missing indicator is used ('separate.missing.indicator == FALSE'), 'na.rm' must be TRUE if missing values occur.\n")}
                              datL <- datL[-isNa,]
                         }
                    }                                                           
                    frml<- as.formula ( paste("~ ",chk[["dependent"]]," - 1",sep="") )
                    datL[, chk[["dependent"]] ] <- as.character( datL[, chk[["dependent"]] ] )
                    matr<- data.frame ( model.matrix ( frml, data = datL) )     
                    datL<- data.frame ( datL,  matr)                            
                    ret <- lapply ( colnames(matr), FUN = function ( dpd ) {
                           attr(datL, "modus") <- modus
                           attr(datL,"depOri") <- chk[["dependent"]]
                           attr(datL,"fc") <- "repTable"
                           res <- repMean ( datL = datL, ID=chk[["ID"]], wgt=chk[["wgt"]], type = type, PSU = chk[["PSU"]], repInd = chk[["repInd"]], repWgt = repWgt,
                                             nest = chk[["nest"]], imp = chk[["imp"]], groups = groups, group.splits = group.splits, group.differences.by=group.differences.by, cross.differences=cross.differences,
                                             crossDiffSE = crossDiffSE, nBoot = nBoot, group.delimiter =group.delimiter, trend = chk[["trend"]], linkErr = chk[["linkErr"]], dependent = dpd, na.rm=na.rm,
                                             doCheck = doCheck, engine=engine, scale = scale, rscales = rscales, mse=mse, rho=rho, verbose=verbose, progress=progress, clusters = NULL)
                           return(res) } )
                    if ( length(ret[[1]][["resT"]]) == 1) {
                         lst <- do.call("rbind", lapply(ret, FUN = function ( x ) { x[["resT"]][["noTrend"]]}))
                         ret <- ret[[1]]
                         ret[["resT"]][["noTrend"]] <- lst
                    }  else  {
                         nams<- names(ret[[1]][["resT"]])
                         lst <- lapply(nams, FUN = function ( na ) { do.call("rbind", lapply(ret, FUN = function ( x ) { x[["resT"]][[na]]}))})
                         names(lst) <- nams
                         ret <- ret[[1]]
                         ret[["resT"]] <- lst
                    }
                    return(ret)
                 }  else  {
                    ret <- eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "table",
                           nest = nest, imp = imp, groups = groups, group.splits = group.splits, group.differences.by = group.differences.by, cross.differences=cross.differences, correct = correct,
                           trend = trend, linkErr = linkErr, dependent = dependent, group.delimiter=group.delimiter, separate.missing.indicator=separate.missing.indicator,
                           expected.values=expected.values, na.rm=na.rm, doCheck=doCheck, modus=modus, engine=engine, scale = scale, rscales = rscales, mse=mse, rho=rho,
                           verbose=verbose, progress=progress, clusters = NULL, fc="repTable")
                 }
                 return(ret)}}


repQuantile<- function(datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"),
            PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL, group.splits = length(groups),
            cross.differences = FALSE, group.delimiter = "_", trend = NULL, linkErr = NULL, dependent, probs = c(0.25, 0.50, 0.75),  na.rm = FALSE,
            nBoot = NULL, bootMethod = c("wSampling","wQuantiles") , doCheck = TRUE, engine = c("survey", "BIFIEsurvey"), 
            scale = 1, rscales = 1, mse=TRUE, rho=NULL, verbose = TRUE, progress = TRUE)  {
            modus      <- identifyMode ( name = "quantile", type = car::recode(match.arg(arg = toupper(type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'"))
            bootMethod <- match.arg ( bootMethod )
            datL       <- eatTools::makeDataFrame ( datL)
            eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "quantile",
                   nest = nest, imp = imp, groups = groups, group.splits = group.splits, cross.differences=cross.differences, trend = trend, linkErr = linkErr, dependent = dependent,
                   group.delimiter=group.delimiter, probs=probs, na.rm=na.rm, nBoot=nBoot, bootMethod=bootMethod, doCheck=doCheck, modus=modus, engine=engine, scale = scale, rscales = rscales, mse=mse, rho=rho, verbose=verbose, progress=progress, clusters=NULL)}


repGlm  <- function(datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"),
            PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL,
            group.splits = length(groups), group.delimiter = "_", cross.differences = FALSE, trend = NULL, linkErr = NULL, formula,
            family=gaussian, forceSingularityTreatment = FALSE, glmTransformation = c("none", "sdY"), doCheck = TRUE, na.rm = FALSE,
            poolMethod = c("mice", "scalar") , useWec = FALSE, scale = 1, rscales = 1, mse=TRUE, rho=NULL,
            hetero=TRUE, se_type = c("HC3", "HC0", "HC1", "HC2", "CR0", "CR2"), clusters = NULL, crossDiffSE.engine= c("lavaan", "lm"), stochasticGroupSizes = FALSE, verbose = TRUE,
            progress = TRUE) {
            datL   <- eatTools::makeDataFrame ( datL)
            modus  <- identifyMode ( name = "glm", type = car::recode(match.arg(arg = toupper(type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'") )
            poolMethod <- match.arg(poolMethod)
            crossDiffSE.engine <- match.arg(crossDiffSE.engine)
            se_type <- chooseSeType(se_type, clusters)
            eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "glm",
                   nest = nest, imp = imp, groups = groups, group.splits = group.splits, cross.differences = cross.differences, trend = trend, linkErr = linkErr,
                   formula=formula, family=family, forceSingularityTreatment=forceSingularityTreatment, glmTransformation = glmTransformation,
                   group.delimiter=group.delimiter, na.rm=na.rm, doCheck=doCheck, modus=modus, poolMethod=poolMethod, useWec=useWec,
                   scale = scale, rscales = rscales, mse=mse, rho=rho, hetero=hetero, se_type=se_type, crossDiffSE.engine=crossDiffSE.engine,
                   stochasticGroupSizes=stochasticGroupSizes, verbose=verbose, progress=progress, clusters=clusters, engine="survey")}

repLmer  <- function(datL, ID, wgt = NULL, L1wgt=NULL, L2wgt=NULL, type = c("JK2", "JK1"),
            PSU = NULL, repInd = NULL, imp=NULL, group=NULL, trend = NULL,  dependent, formula.fixed, formula.random,
            doCheck = TRUE, na.rm = FALSE, clusters, verbose = TRUE) {
            datL   <- eatTools::makeDataFrame ( datL)
            modus  <- identifyMode ( name = "lmer", type = car::recode(match.arg(arg = toupper(type), choices = c("JK2", "JK1")), "'FAY'='Fay'") )
            eatRep(datL =datL, ID=ID , wgt = wgt, L1wgt=L1wgt, L2wgt=L2wgt, type=type, PSU = PSU, repInd = repInd, toCall = "lmer",
                   imp = imp, groups =group, trend = trend, dependent=dependent, formula.fixed=formula.fixed, formula.random=formula.random,engine="BIFIEsurvey",
                   na.rm=na.rm, doCheck=doCheck, modus=modus, verbose=verbose, clusters=clusters)}

eatRep <- function (datL, ID, wgt = NULL, L1wgt=NULL, L2wgt=NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"), PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL,
          toCall = c("mean", "table", "quantile", "glm", "cov", "lmer", "glmer"), groups = NULL, refGrp = NULL, group.splits = length(groups), group.differences.by = NULL,
          cross.differences = FALSE, group.delimiter = "_", adjust=NULL, useEffectLiteR = TRUE, trend = NULL, linkErr = NULL, dependent, na.rm = FALSE, forcePooling = TRUE, boundary = 3, doCheck = TRUE,
          separate.missing.indicator = FALSE, expected.values = NULL, probs = NULL, nBoot = NULL, bootMethod = NULL, formula=NULL, family=NULL, formula.fixed=NULL, formula.random=NULL,
          forceSingularityTreatment = FALSE, glmTransformation = c("none", "sdY"), correct, onlyCheck = FALSE, modus, poolMethod = "mice", useWec = FALSE, engine,
          scale, rscales, mse, rho, reihenfolge = NULL, hetero, se_type, crossDiffSE.engine, stochasticGroupSizes, verbose, progress, clusters=NULL, fc = NULL, isRecursive = FALSE, depOri = NULL) {
          datL  <- eatTools::makeDataFrame(datL, name = "datL")
          if ( isTRUE(useWec) ) { forceSingularityTreatment <- TRUE; poolMethod <- "scalar"}
          if (is.null(fc) && isFALSE(onlyCheck)) {                              
#               beg   <- Sys.time()
               fc    <- identifyFunctionCall()
#               message(paste0("Identify function call: ", timeFormat(Sys.time() - beg)))
          }
          toCall<- match.arg(toCall)                                            
          type  <- car::recode(match.arg(arg = toupper(type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'")
          if ( type == "NONE") {doJK <- FALSE }  else {doJK <- TRUE}
          engine<- match.arg(arg = engine, choices = c("survey", "BIFIEsurvey"))
          glmTransformation <- match.arg(glmTransformation)
          if(isFALSE(forceSingularityTreatment) & glmTransformation != "none") {
             message("'forceSingularityTreatment' was set to 'FALSE'. Please note that 'glmTransformation' is only possible if 'forceSingularityTreatment' is 'TRUE'.")
          }
          uv_av <- identify_UV_AV(toCall=toCall, dependent=dependent, formula=formula, formula.random=formula.random, formula.fixed=formula.fixed, glmerFormula=NULL)
          if(is.null(groups))  {                                                
             groups <- "wholeGroup"
             datL[,"wholeGroup"] <- "wholePop"
             groupWasNULL <- TRUE
          }  else  {
             groupWasNULL <- FALSE
          }
          if ( length(linkErr) > 1 ) {
             leFrame <- linkErr
             linkErr <- NULL
          }  else  {
             leFrame <- NULL
          }
          allVar<- list(ID = ID, wgt = wgt, L1wgt=L1wgt, L2wgt=L2wgt, wgt, PSU = PSU, repInd = repInd, repWgt = repWgt, nest=nest, imp=imp, group = groups, trend=trend, linkErr = linkErr, group.differences.by=group.differences.by, dependent = uv_av[["av"]], independent=uv_av[["uv"]], adjust=adjust, clusters=clusters)
          allNam<- lapply(allVar, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = datL, variable=ii, warnIfMissing = TRUE)})
          if ( !is.null(leFrame)){
             allNam[["linkErr"]] <- leFrame
          }
          if(forceSingularityTreatment == TRUE && !is.null(allNam[["PSU"]]) ) { poolMethod <- "scalar"}
          foo   <- checkJK.arguments(type=type, repWgt=repWgt, PSU=PSU, repInd=repInd)
          auchUV<- checkWecForUV(dat=datL, allNam = allNam)
          if (isFALSE(isRecursive)) {                                           
#              beg   <- Sys.time()
              datL  <- checkGroupVars ( datL = datL, allNam = allNam, auchUV = auchUV)
#              message(paste0("checkGroupVars: ", timeFormat(Sys.time() - beg)))
          }
          datNam<- checkForAdjustmentAndLmer (datL=datL, allNam=allNam, groupWasNULL=groupWasNULL, formula.random=formula.random, formula.fixed=formula.fixed)
          formula.fixed <- datNam[["formula.fixed"]]                            
          formula.random <- datNam[["formula.random"]]                          
          foo   <- checkNameConvention( allNam = datNam[["allNam"]])
          if (isTRUE(useWec) ) {
              if ( length(allNam[["independent"]]) != 1 ) {stop("Only one independent (grouping) variable is allowed for weighted effect coding.\n")}
              if(!inherits(datNam[["datL"]][,allNam[["independent"]]],c("factor", "character", "logical", "integer"))) {stop(paste0("For weighted effect coding, independent (grouping) variable '",allNam[["independent"]],"' must be of class 'factor', 'character', 'logical', or 'integer'.\n"))}
          }
          if (poolMethod == "mice" && !is.null(datNam[["allNam"]][["nest"]]) && toCall == "glm" ) {
              message("Method 'mice' is not available for nested imputation. Switch to method 'scalar'.")
              poolMethod <- "scalar"
          }
          engine<- checkEngine (engine = engine, allNam=datNam[["allNam"]], modus=modus, toCall = toCall, type=type)
          if ( isTRUE(onlyCheck) ) {
              ret <- datNam[["allNam"]]
          }  else  {
              if(!is.null(datNam[["allNam"]][["trend"]])) {                     
                  nlev<- length(unique(datNam[["datL"]][,datNam[["allNam"]][["trend"]]]))
                  if(nlev == 1) {stop(paste0("Trend variable '",datNam[["allNam"]][["trend"]],"' is constant with value '",unique(datNam[["datL"]][,datNam[["allNam"]][["trend"]]]),"'."))}
                  chk5<- checkFactorLevels(datL=datNam[["datL"]], allNam=datNam[["allNam"]])
                  resT<- by ( data = datNam[["datL"]], INDICES = datNam[["datL"]][,datNam[["allNam"]][["trend"]]], FUN = function ( subdat ) {
                         if(verbose) {cat(paste("\nTrend group: '",subdat[1,datNam[["allNam"]][["trend"]] ], "'\n",sep=""))}
                         do    <- paste ( "eatRep ( ", paste(names(formals(eatRep)), car::recode(names(formals(eatRep)), "'trend'='NULL'; 'datL'='subdat'; 'group.differences.by'='group.differences.by'; 'isRecursive'='TRUE'"), sep =" = ", collapse = ", "), ")",sep="")
                         foo   <- eval(parse(text=do))
                         foo[["resT"]][["noTrend"]][,datNam[["allNam"]][["trend"]]] <- subdat[1,datNam[["allNam"]][["trend"]] ]
                         out1  <- foo[["resT"]][["noTrend"]]
                         out2  <- foo[["allNam"]]
                         return(list(out1=out1, out2=out2))}, simplify = FALSE)
                  allNam <- c(datNam[["allNam"]], resT[[1]][["out2"]])
                  allNam <- allNam[unique(names(allNam))]
                  le     <- createLinkingError ( allNam = allNam, resT = resT, datL = datNam[["datL"]], fc=fc, toCall = toCall)
                  out3   <- lapply(resT, FUN = function ( k ) { k[["out1"]]})
                  ret    <- list(resT = out3, allNam = allNam, toCall = toCall, family=family, le=le)
                  return(ret)
              }  else {
                  if( length( setdiff ( datNam[["allNam"]][["group.differences.by"]],datNam[["allNam"]][["group"]])) != 0) {stop("Variable in 'group.differences.by' must be included in 'groups'.\n")}
                  toAppl<- superSplitter(group = datNam[["allNam"]][["group"]], group.splits = group.splits, group.differences.by = datNam[["allNam"]][["group.differences.by"]], group.delimiter = group.delimiter , dependent=datNam[["allNam"]][["dependent"]] )
                  if(verbose){cat(paste(length(toAppl)," analyse(s) overall according to: 'group.splits = ",paste(group.splits, collapse = " ") ,"'.", sep=""))}
                  ret   <- createAnalysisInfTable(toAppl=toAppl, verbose=verbose, allNam=datNam[["allNam"]])
                  allNam<- setCrossDifferences (cross.differences=cross.differences, allNam=datNam[["allNam"]], group.splits=group.splits)
#                  beg   <- Sys.time()
                  cls   <- createLoopStructure(datL = datNam[["datL"]], allNam = allNam, verbose=verbose)
#                  message(paste0("createLoopStructure: ", timeFormat(Sys.time() - beg)))
                  if(!is.null(cls[["allNam"]][["cross.differences"]])) {
                      if(length(cls[["allNam"]][["group"]])>1) {
                         lev <- unlist(lapply(cls[["allNam"]][["group"]], FUN = function ( v ) { unique(as.character(cls[["datL"]][,v]))}))
                         if (length(lev) != length(unique(lev))) {stop("Factor levels of grouping variables are not disjunct.\n")}
                      }
                  }
                  if(toCall %in% c("mean", "quantile", "glm")) {
                     if(!inherits(cls[["datL"]][,cls[["allNam"]][["dependent"]]],  c("integer", "numeric"))) {
                         stop(paste0("Dependent variable '",cls[["allNam"]][["dependent"]],"' has to be of class 'integer' or 'numeric'.\n"))
                     }
                  }
                  repA  <- assignReplicates ( repWgt=repWgt, allNam=cls[["allNam"]], datL = cls[["datL"]], engine = engine, type=type , progress=progress, verbose=verbose)
                  allRes<- innerLoop(toAppl=toAppl, ret=ret, toCall=toCall, allNam=cls[["allNam"]], datL=cls[["datL"]], separate.missing.indicator=separate.missing.indicator,
                           na.rm=na.rm, expected.values=expected.values, group.delimiter=group.delimiter,type=type, repA=repA, modus=modus,probs=probs, nBoot=nBoot,
                           bootMethod=bootMethod, formula=formula, forceSingularityTreatment=forceSingularityTreatment,glmTransformation=glmTransformation,doJK=doJK,
                           poolMethod=poolMethod, useWec=useWec,refGrp=refGrp, scale = scale,rscales = rscales, mse=mse, rho=rho,reihenfolge=reihenfolge, hetero=hetero, se_type=se_type,
                           useEffectLiteR=useEffectLiteR,crossDiffSE.engine=crossDiffSE.engine,stochasticGroupSizes=stochasticGroupSizes, progress=progress, correct=correct,family=family,
                           verbose=verbose, doCheck=doCheck, engine=engine, L1wgt=L1wgt, L2wgt=L2wgt,formula.fixed=formula.fixed, formula.random=formula.random )
                  if(verbose){cat("\n")}
                  allRes <- clearTab(allRes, allNam = cls[["allNam"]], depVarOri = depOri, fc=fc, toCall=toCall, datL = cls[["datL"]])
                  allRes <- list(resT = list(noTrend = allRes), allNam = cls[["allNam"]], toCall = toCall, family=family)
                  return(allRes) }} }
                  
checkFactorLevels <- function(datL, allNam) {
       if (!is.null(allNam[["group"]])) {
             foo <- lapply(allNam[["group"]], FUN = function ( gr ) {
                    ch <- by(data = datL, INDICES = datL[,allNam[["trend"]]], FUN = function ( subdat ) { table(subdat[,gr]) }, simplify = FALSE )
                    cmb<- combinat::combn(x=1:length(ch), m=2, simplify=FALSE)
                    ch1<- all(unlist(lapply(cmb,FUN = function (y) {all ( names(ch[[y[1]]]) == names(ch[[y[2]]]))})))
                    if(!ch1) {
                        tab <- table(datL[,c(gr,allNam[["trend"]])])
                        stop(paste0("Levels of grouping variable '",gr,"' do not match between trend groups: \n", eatTools::print_and_capture(tab, 5) ))
                    }})}}

createAnalysisInfTable <- function(toAppl, verbose, allNam) {
         if ( length ( toAppl ) > 1) {
               ret <- do.call("rbind", lapply(1:length(toAppl), FUN = function ( y ) {
               gdb <- attr(toAppl[[y]], "group.differences.by")
               if ( is.null(gdb)) {gdb <- NA}
               res <-  data.frame ( analysis.number = y, hierarchy.level = length(toAppl[[y]]), groups.divided.by = paste(toAppl[[y]], collapse=" + "), group.differences.by = gdb)
               if ( !is.null(allNam[["adjust"]])) {
                      res[,"adjust"] <- car::recode(res[,"hierarchy.level"], "0='FALSE'; else = 'TRUE'")
               }
               return(res)}))
               if(verbose){cat("\n \n"); print(ret, row.names=FALSE)}
         }  else  {ret <- NULL}
         return(ret)}

innerLoop <- function (toAppl, ret, toCall, allNam, datL, separate.missing.indicator, na.rm, expected.values, group.delimiter,type, repA, modus,probs,nBoot,bootMethod, formula, forceSingularityTreatment,
                      glmTransformation,doJK, poolMethod, useWec,refGrp, scale,rscales, mse, rho,reihenfolge, hetero, se_type, useEffectLiteR,crossDiffSE.engine,stochasticGroupSizes, progress, correct,
                      family, verbose, doCheck, engine, L1wgt, L2wgt,formula.fixed, formula.random)  {
       allRes<- do.call(plyr::rbind.fill, lapply( names(toAppl), FUN = function ( gr ) {
                str1  <- createInfoString (ai = ret, toCall=toCall, gr=gr, toAppl=toAppl)
                if(toCall %in% c("mean", "table"))  { allNam[["group.differences.by"]] <- attr(toAppl[[gr]], "group.differences.by") }
                if( nchar(gr) == 0 ){                                           
                    datL[,"dummyGroup"] <- "wholeGroup"
                    allNam[["group"]] <- "dummyGroup"                           ### problematisch!! "allNam" wird hier in jedem Schleifendurchlauf ueberschrieben -- nicht so superschoen!
                    allNam[["adjust"]] <- NULL                                  
                } else {
                    allNam[["group"]] <- toAppl[[gr]]
                }
                noMis <- unlist ( c ( allNam[-na.omit(match(c("group", "dependent", "cross.differences"), names(allNam)))], toAppl[gr]) )
                miss  <- which ( sapply(datL[,noMis], FUN = function (uu) {length(which(is.na(uu)))}) > 0 )
                if(length(miss)>0) { warning("Unexpected missings in variable(s) ",paste(names(miss), collapse=", "),".")}
#                beg   <- Sys.time()
                datL  <- checkImpNest(datL = datL, doCheck=doCheck, toAppl = toAppl, gr=gr, allNam = allNam, toCall=toCall, separate.missing.indicator=separate.missing.indicator, na.rm=na.rm)
#                message(paste0("checkImpNest: ", timeFormat(Sys.time() - beg)))
                pev   <- prepExpecVal (toCall = toCall, expected.values=expected.values, separate.missing.indicator=separate.missing.indicator, allNam=allNam, datL = datL)
                if ( engine=="survey" || isFALSE(doJK)) {
                anaA<- do.call("rbind", by(data = pev[["datL"]], INDICES = pev[["datL"]][,"isClear"], FUN = doSurveyAnalyses, allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter,
                       type=type, repA=repA, modus=modus, separate.missing.indicator = separate.missing.indicator, expected.values=pev[["expected.values"]], probs=probs,
                       nBoot=nBoot,bootMethod=bootMethod, formula=formula, forceSingularityTreatment=forceSingularityTreatment, glmTransformation=glmTransformation,
                       toCall=toCall, doJK=doJK, poolMethod=poolMethod, useWec=useWec, refGrp=refGrp, scale = scale, rscales = rscales, mse=mse, rho=rho,
                       reihenfolge=reihenfolge, hetero=hetero, se_type=se_type, useEffectLiteR=useEffectLiteR, crossDiffSE.engine=crossDiffSE.engine,
                       stochasticGroupSizes=stochasticGroupSizes, progress=progress, correct=correct, family=family, str1=str1))
                }  else  {                                                      
                       anaA<- do.call("rbind", by(data = pev[["datL"]], INDICES = pev[["datL"]][,"isClear"], FUN = doBifieAnalyses, allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter,
                              separate.missing.indicator = separate.missing.indicator, expected.values=pev[["expected.values"]], probs=probs, formula=formula, glmTransformation=glmTransformation,
                              toCall=toCall, modus=modus, type=type, verbose=verbose, L1wgt=L1wgt, L2wgt=L2wgt,formula.fixed=formula.fixed, formula.random=formula.random))
                }
                if( "dummyGroup" %in% colnames(anaA) )  { anaA <- anaA[,-match("dummyGroup", colnames(anaA))] }
                return(anaA)}))
       rownames(allRes) <- NULL
       return(allRes)}


identifyFunctionCall <- function() {
          i     <- 0                                                            
          fc    <- NULL                                                         
          while ( !is.null(sys.call(i))) { fc <- c(fc, eatTools::crop(unlist(strsplit(deparse(sys.call(i))[1], split = "\\("))[1])); i <- i-1  }
          fc   <- unlist(lapply(strsplit(fc, ":"), FUN = function ( l ) {l[length(l)]}))
          fc   <- fc[max(which(fc %in% c("repMean", "repTable", "repGlm", "repQuantile", "repLmer", "repGlmer")))]
          return(fc)}

checkRegression <- function ( dat, allNam, useWec ) {
                   ch <- lapply( allNam[["independent"]], FUN = function ( i ) {
                         isKonst <- length(unique(dat[,i]))
                         if ( isKonst == 1) {
                              warning("Predictor '",i,"' is constant. Please check your data.")
                         }
                         if ( inherits ( dat[,i],  "character" )) {
                              warning("Predictor '",i,"' has class 'character'. Please check your data.")
                         }
                         if ( inherits ( dat[,i],  c("character", "factor") )) {
                              if ( isKonst > 15 && isFALSE(useWec) ) {
                                   warning("Predictor '",i,"' of class '",class ( dat[,i] ),"' has ",isKonst," levels. Please check whether this is intended.")
                              }
                         } })   }                                               

createLinkingError <- function  ( allNam = allNam, resT = resT, datL = datL, fc, toCall) {
          if (length(allNam[["linkErr"]]) > 1) {                                
              le <- allNam[["linkErr"]]
              attr(le, "linkingErrorFrame") <- TRUE                             
              return(le)
          }  else  {
              if ( is.null ( allNam[["linkErr"]] ) ) {
                   message("Note: No linking error was defined. Linking error will be defaulted to '0'.")
                   allNam[["linkErr"]] <- "le"
                   datL[,"le"]         <- 0
              }
              init <- data.frame ( trendVar = allNam[["trend"]], trendLevel1 = sort(unique(datL[,allNam[["trend"]]]))[1] , trendLevel2 = sort(unique(datL[,allNam[["trend"]]]))[2], stringsAsFactors = FALSE)
              if ( fc == "repTable") {
                   le <- data.frame ( init, depVar = unique(resT[[1]][["out1"]][,"depVar"]), unique(datL[, c(unique(resT[[1]][["out1"]][,"depVar"]), allNam[["linkErr"]]),drop=FALSE]), stringsAsFactors = FALSE)
              }
              if ( fc == "repMean") {
                   stopifnot (length(unique(datL[,allNam[["linkErr"]]])) == 1)
                   le <- unique(datL[,allNam[["linkErr"]]])                     
                   le <- data.frame ( init, depVar = unique(resT[[1]][["out1"]][,"depVar"]), parameter = c("mean", "sd"), le = c(le,0), stringsAsFactors = FALSE)
              }
              if ( fc %in% c("repGlm", "repLmer", "repGlmer")) {
                   stopifnot (length(unique(datL[,allNam[["linkErr"]]])) == 1)
                   le <- unique(datL[,allNam[["linkErr"]]])
                   le <- data.frame ( init, depVar = unique(as.character(resT[[1]][["out1"]][,"depVar"])), parameter = unique(as.character(resT[[1]][["out1"]][,"parameter"])), le = le, stringsAsFactors = FALSE)
              }
              if ( fc == "repQuantile") {
                   stopifnot (length(unique(datL[,allNam[["linkErr"]]])) == 1)
                   le <- data.frame ( init, unique(resT[[1]][["out1"]][,c("depVar", "parameter")]), le = unique(datL[,allNam[["linkErr"]]]), stringsAsFactors = FALSE)
              }
              colnames(le)[5:6] <- c("parameter", "le")
              if ( nrow(le) > length(unique(le[,"parameter"]))) {
                   stop("Linking errors must be unique for levels of dependent variable.\n")
              }
          }
          return(le)}

conv.quantile      <- function ( dat.i , allNam, na.rm, group.delimiter, probs, nBoot,bootMethod, modus) {
                      ret  <- do.call("rbind", by(data = dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( sub.dat) {
                              if( all(sub.dat[,allNam[["wgt"]]] == 1) )  {      ### alle Gewichte sind 1 bzw. gleich
                                 ret   <- Hmisc::hdquantile(x = sub.dat[,allNam[["dependent"]]], se = TRUE, probs = probs,na.rm=na.rm )
                                 ret   <- data.frame (group = paste(sub.dat[1,allNam[["group"]],drop=FALSE], collapse=group.delimiter), depVar = allNam[["dependent"]], modus = modus, parameter = rep(names(ret),2), coefficient = rep(c("est","se"),each=length(ret)),value = c(ret,attr(ret,"se")),sub.dat[1,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                              } else {                                          ### wenn Gewichte gefordert, koennen SEs ueber Bootstrap bestimmt werden
                                 if(!is.null(nBoot)) {
                                     if(nBoot<5) {nBoot <- 5}
                                     if(bootMethod == "wQuantiles") {           ### Variante 1
                                         x     <- sub.dat[,allNam[["dependent"]]]
                                         ret   <- boot::boot(data = x, statistic = function ( x, i) {Hmisc::wtd.quantile(x = x[i], weights = sub.dat[i,allNam[["wgt"]]], probs = probs,na.rm=na.rm )}, R=nBoot)
                                         ret   <- data.frame (group = paste(sub.dat[1,allNam[["group"]],drop=FALSE], collapse=group.delimiter), depVar = allNam[["dependent"]], modus = modus, parameter = rep(as.character(probs),2), coefficient = rep(c("est","se"),each=length(probs)), value = c(ret$t0, sapply(data.frame(ret$t), sd)), sub.dat[1,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                                     } else {                                   ### Variante 2
                                         ret   <- do.call("rbind", lapply(1:nBoot, FUN = function (b){
                                                  y   <- sample(x = sub.dat[,allNam[["dependent"]]], size = length(sub.dat[,allNam[["dependent"]]]), replace = TRUE, prob = sub.dat[,allNam[["wgt"]]]/sum(sub.dat[,allNam[["wgt"]]]))
                                                  ret <- Hmisc::hdquantile(x = y, se = FALSE, probs = probs,na.rm=na.rm )
                                                  return(ret)}))
                                         ret   <- data.frame (group = paste(sub.dat[1,allNam[["group"]],drop=FALSE], collapse=group.delimiter), depVar = allNam[["dependent"]], modus = modus, parameter = rep(as.character(probs),2), coefficient = rep(c("est","se"),each=length(probs)), value = c(Hmisc::wtd.quantile(x = sub.dat[,allNam[["dependent"]]], weights = sub.dat[,allNam[["wgt"]]], probs = probs,na.rm=na.rm ), sapply(data.frame(ret),sd)) , sub.dat[1,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                                     }
                                 } else {
                                     ret   <- Hmisc::wtd.quantile(x = sub.dat[,allNam[["dependent"]]], weights = sub.dat[,allNam[["wgt"]]], probs = probs,na.rm=na.rm )
                                     ret   <- data.frame (group = paste(sub.dat[1,allNam[["group"]],drop=FALSE], collapse=group.delimiter), depVar = allNam[["dependent"]], modus = modus, parameter = rep(as.character(probs),2), coefficient = rep(c("est","se"),each=length(probs)), value = c(ret, rep(NA, length(probs))) , sub.dat[1,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                                 }
                              }
                              return(ret)}))
                      ret[,"comparison"] <- NA
                      return(eatTools::facToChar(ret))}

jackknife.quantile <- function ( dat.i , allNam, na.rm, type, repA, probs, group.delimiter, modus, scale , rscales, mse, rho) {
                      typeS          <- car::recode(type, "'JK2'='JKn'")        ### typeS steht fuer type_Survey
                      design         <- svrepdesign(data = dat.i[,c(allNam[["group"]], allNam[["dependent"]]) ], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                      formel         <- as.formula(paste("~ ",allNam[["dependent"]], sep = "") )
    ### Hotfix: return.replicates = FALSE gesetzt, weil es sonst ab survey version 4.1-1 eine fehlermeldung gibt ... weiss nicht, ob man die replicates spaeter nochmal braucht, ich glaube nicht
                      quantile.imp   <- svyby(formula = formel, by = as.formula(paste("~", paste(allNam[["group"]], collapse = " + "))), design = design, FUN = svyquantile, quantiles = probs, return.replicates = FALSE, na.rm = na.rm)
                      molt           <- eatTools::facToChar(reshape2::melt(data=quantile.imp, id.vars=allNam[["group"]], na.rm=FALSE))
                      molt[,"parameter"]  <- eatTools::crop(eatTools::removePattern(eatTools::removePattern(molt[,"variable"],paste0("se.", allNam[["dependent"]])),allNam[["dependent"]]), char=".")
                      molt           <- do.call("rbind", plyr::alply(molt, .margins = 1, .fun = function (zeile) {
                                        coef <- eatTools::removePattern(eatTools::crop(eatTools::removePattern(zeile[["variable"]], zeile[["parameter"]]), "."), allNam[["dependent"]])
                                        if ( coef=="") {coef <- "est"} else {stopifnot(coef == "se."); coef <- "se"}
                                        zeile[["coefficient"]] <- coef
                                        return(zeile)}))
                      return(eatTools::facToChar(data.frame ( group = apply(molt[,allNam[["group"]],drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = allNam[["dependent"]], modus = paste(modus,"survey", sep="__"), comparison = NA, molt[,c("parameter", "coefficient", "value", allNam[["group"]])], stringsAsFactors = FALSE))) }


conv.table      <- function ( dat.i , allNam, na.rm, group.delimiter, separate.missing.indicator , correct, expected.values, modus) {
                   table.cast <- do.call("rbind", by(data = dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( sub.dat) {
                                 prefix <- data.frame(sub.dat[1,allNam[["group"]], drop=FALSE], row.names = NULL, stringsAsFactors = FALSE )
                                 foo    <- make.indikator(variable = sub.dat[,allNam[["dependent"]]], name.var = "ind", force.indicators =expected.values, separate.missing.indikator = "no")
                                 if(all(dat.i[,allNam[["wgt"]]] == 1)) {ret    <- data.frame ( prefix , eatTools::descr(foo[,-1, drop = FALSE],na.rm=TRUE)[,c("Mean", "std.err")], stringsAsFactors = FALSE )
                                 } else { ret    <- data.frame ( prefix , eatTools::descr(foo[,-1, drop = FALSE], p.weights = sub.dat[,allNam[["wgt"]]],na.rm=TRUE)[,c("Mean", "std.err")], stringsAsFactors = FALSE )}
                                 ret[,"parameter"] <- substring(rownames(ret),5)
                                 return(ret)}) )
                   Ns         <- do.call("rbind", by(dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function (y) {data.frame ( y[1,allNam[["group"]], drop=FALSE], parameter = "Ncases", Mean=length(unique(y[,allNam[["ID"]]])), stringsAsFactors = FALSE) }))
                   table.cast <- plyr::rbind.fill(table.cast, Ns)
                   if(!is.null(allNam[["group.differences.by"]]))   {
                      m            <- table.cast
                      m$comb.group <- apply(m, 1, FUN = function (ii) { eatTools::crop(paste( ii[allNam[["group"]]], collapse = "."))})
                      m$all.group  <- 1
                      res.group    <- tempR <- setdiff(allNam[["group"]], allNam[["group.differences.by"]])
                      if(length(res.group) == 0 ) {res.group <- "all.group"}
                      difs         <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                      if(length(tempR)>0) {
                                         datSel    <- merge(dat.i, iii[!duplicated(iii[,res.group]),res.group,drop=FALSE], by = res.group, all = FALSE)
                                      }  else  {
                                         datSel    <- dat.i
                                      }
                                      tbl    <- table(datSel[,c(allNam[["group.differences.by"]], allNam[["dependent"]])])
                                      chisq  <- chisq.test(tbl, correct = correct)
                                      scumm  <- iii[!duplicated(iii[,res.group]),res.group,drop = FALSE]
                                      group  <- paste( paste( colnames(scumm), as.character(scumm[1,]), sep="="), sep="", collapse = ", ")
                                      dif.iii<- data.frame(group = group, parameter = "chiSquareTest", comparison = "groupDiff", depVar = allNam[["dependent"]], modus=modus, coefficient = c("chi2","df","pValue"), value = c(chisq[["statistic"]],chisq[["parameter"]],chisq[["p.value"]]) , stringsAsFactors = FALSE )
                                      return(dif.iii)}))                        
                   }                                                            
                   ret        <- reshape2::melt(table.cast, measure.vars = c("Mean", "std.err"), na.rm=TRUE)
                   ret[,"coefficient"] <- car::recode(ret[,"variable"], "'Mean'='est'; 'std.err'='se'")
                   ret        <- data.frame ( group = apply(ret[,allNam[["group"]],drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = allNam[["dependent"]], modus = modus, comparison = NA, ret[,c("coefficient", "parameter")], value = ret[,"value"], ret[,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                   if(!is.null(allNam[["group.differences.by"]]))   {return(eatTools::facToChar(plyr::rbind.fill(ret,difs)))} else {return(eatTools::facToChar(ret))}}


jackknife.table <- function ( dat.i , allNam, na.rm, group.delimiter, type, repA, separate.missing.indicator, expected.values, modus, scale, rscales, mse, rho) {
                   dat.i[,allNam[["dependent"]]] <- factor(dat.i[,allNam[["dependent"]]], levels = expected.values)
                   typeS     <- car::recode(type, "'JK2'='JKn'")
                   design    <- svrepdesign(data = dat.i[,c(allNam[["group"]], allNam[["dependent"]])], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                   formel    <- as.formula(paste("~factor(",allNam[["dependent"]],", levels = expected.values)",sep=""))
                   means     <- svyby(formula = formel, by = as.formula(paste("~", paste(as.character(allNam[["group"]]), collapse = " + "))), design = design, FUN = svymean, deff = FALSE, return.replicates = TRUE)
                   Ns        <- do.call("rbind", by(dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function (y) {data.frame ( y[1,allNam[["group"]], drop=FALSE], variable = "est____________Ncases", value=length(unique(y[,allNam[["ID"]]])), stringsAsFactors = FALSE) }))
                   cols      <- match(paste("factor(",allNam[["dependent"]],", levels = expected.values)",expected.values,sep=""), colnames(means))
                   colnames(means)[cols] <- paste("est",expected.values, sep="____________")
                   cols.se   <- grep("^se[[:digit:]]{1,5}$", colnames(means) )
                   stopifnot(length(cols) == length(cols.se))
                   colnames(means)[cols.se] <- paste("se____________", expected.values, sep="")
                   molt      <- reshape2::melt(data=means, id.vars=allNam[["group"]], na.rm=TRUE)
                   molt      <- rbind(molt, Ns)
                   splits    <- data.frame ( do.call("rbind", strsplit(as.character(molt[,"variable"]),"____________")), stringsAsFactors = FALSE)
                   colnames(splits) <- c("coefficient", "parameter")
                   ret       <- data.frame ( group = apply(molt[,allNam[["group"]],drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = allNam[["dependent"]], modus = paste(modus,"survey", sep="__"), comparison = NA, splits, value = molt[,"value"], molt[,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                   if(!is.null(allNam[["group.differences.by"]]))   {
                      m            <- ret
                      m$comb.group <- apply(m, 1, FUN = function (ii) { eatTools::crop(paste( ii[allNam[["group"]]], collapse = "."))})
                      m$all.group  <- 1
                      res.group    <- tempR <- setdiff(allNam[["group"]], allNam[["group.differences.by"]])
                      if(length(res.group) == 0 ) {res.group <- "all.group"}
                      difs           <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                        if(length(tempR)>0) {
                                           datSel    <- merge(dat.i, iii[!duplicated(iii[,res.group]),res.group,drop=FALSE], by = res.group, all = FALSE)
                                        }  else  {
                                           datSel    <- dat.i
                                        }
                                        designSel <- svrepdesign(data = datSel[,c(allNam[["group.differences.by"]], allNam[["dependent"]])], weights = datSel[,allNam[["wgt"]]], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(datSel[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                                        formel    <- as.formula(paste("~", allNam[["group.differences.by"]] ,"+",allNam[["dependent"]],sep=""))
                                        tbl       <- svychisq(formula = formel, design = designSel, statistic = "Chisq")
                                        scumm     <- iii[!duplicated(iii[,res.group]),res.group,drop = FALSE]
                                        group     <- paste( paste( colnames(scumm), as.character(scumm[1,]), sep="="), sep="", collapse = ", ")
                                        dif.iii   <- data.frame(group = group, parameter = "chiSquareTest", depVar = allNam[["dependent"]], modus = paste(modus,"survey", sep="__"), coefficient = c("chi2","df","pValue"), value = c(tbl[["statistic"]],tbl[["parameter"]],tbl[["p.value"]]) , stringsAsFactors = FALSE )
                                        return(dif.iii)                         
                      } ))                                                      
                      difs[,"comparison"] <- "groupDiff"
                   }
                   if(!is.null(allNam[["group.differences.by"]]))   {return(eatTools::facToChar(plyr::rbind.fill(ret,difs)))} else {return(eatTools::facToChar(ret))}}


conv.mean      <- function (dat.i , allNam, na.rm, group.delimiter, modus) {
                  deskr    <- data.frame ( do.call("rbind", by(data = dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( sub.dat) {
                              prefix <- sub.dat[1,allNam[["group"]], drop=FALSE]
                              if ( all(sub.dat[,allNam[["wgt"]]] == 1) )  {
                                   useWGT <- NULL
                              }  else  {
                                   useWGT <- sub.dat[,allNam[["wgt"]]]
                                   attr(useWGT, "onlyUnweightedN") <- TRUE      
                              }                                                 
                              ret    <- data.frame ( nValidUnweighted = length(na.omit(sub.dat[, allNam[["dependent"]] ])), prefix, eatTools::descr(sub.dat[, allNam[["dependent"]] ], p.weights = useWGT, na.rm=na.rm)[,c("N", "N.valid", "Mean", "std.err", "Var", "SD")], stringsAsFactors = FALSE)
                              names(ret) <- c( "nValidUnweighted", allNam[["group"]] , "Ncases", "NcasesValid", "mean", "se.mean", "var","sd")
                              return(ret)})), modus=modus, stringsAsFactors = FALSE)
                  if(!is.null(allNam[["group.differences.by"]]))   {
                     nCat <- table(as.character(dat.i[,allNam[["group.differences.by"]]]))
                     if ( length(nCat) < 2 ) {
                          cat(paste("Warning: Grouping variable '", allNam[["group.differences.by"]], "' only has one category within imputation and/or nest. Group differences cannot be computed. Skip computation.\n",sep=""))
                     }  else  {
                          m            <- deskr
                          m$comb.group <- apply(m, 1, FUN = function (ii) { eatTools::crop(paste( ii[allNam[["group"]]], collapse = "."))})
                          m$all.group  <- 1
                          res.group    <- tempR <- setdiff(allNam[["group"]], allNam[["group.differences.by"]])
                          if(length(res.group) == 0 ) {res.group <- "all.group"}
                          kontraste    <- combinat::combn ( x = sort(unique(as.character(m[,allNam[["group.differences.by"]]]))), m = 2, simplify = FALSE)
                          difs         <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                          ret <- do.call("rbind", lapply(kontraste, FUN = function ( k ) {
                                                 if ( sum ( k %in% iii[,allNam[["group.differences.by"]]]) != length(k) ) {
                                                    warning("Cannot compute contrasts for 'group.differences.by = ",allNam[["group.differences.by"]],"'.")
                                                    return(NULL)
                                                 }  else  {
                                                    vgl.iii   <- iii[iii[,allNam[["group.differences.by"]]] %in% k ,]
                                                    stopifnot ( nrow(vgl.iii) == 2)
                                                    true.diff <- diff(vgl.iii[,"mean"])
                                                    scumm     <- sapply(vgl.iii[,res.group,drop = FALSE], as.character)
                                                    group     <- paste( paste( colnames(scumm), scumm[1,], sep="="), sep="", collapse = ", ")
                                                    dummy     <- do.call("cbind", lapply ( allNam[["group"]], FUN = function ( gg ) {
                                                                 ret <- data.frame ( paste ( unique(vgl.iii[,gg]), collapse = ".vs."))
                                                                 colnames(ret) <- gg
                                                                 return(ret)}))
                                                    dif.iii   <- data.frame(dummy, group = paste(group, paste(k, collapse = ".vs."),sep="____"), parameter = "mean", coefficient = c("est","se"), depVar = allNam[["dependent"]], modus=modus, value = c(true.diff, sqrt( sum(vgl.iii[,"sd"]^2 / vgl.iii[,"nValidUnweighted"]) )) , stringsAsFactors = FALSE )
                                                    stopifnot(nrow(dif.iii)==2, nrow(vgl.iii) == 2)
                                                    dummy2    <- dif.iii[1,]
                                                    dummy2[,"coefficient"] <- "es"
                                                    dummy2[,"value"]       <- dif.iii[which(dif.iii[,"coefficient"] == "est"),"value"] / sqrt(0.5*sum(vgl.iii[,"sd"]^2))
                                                    return(dif.iii)             
                                                 } }))                          
                                          return(ret)}))
                        difs[,"comparison"] <- "groupDiff"
                     }
                  }
                  deskrR   <- reshape2::melt(data = deskr, id.vars = allNam[["group"]], measure.vars = setdiff(colnames(deskr), c("nValidUnweighted", "modus", allNam[["group"]]) ), na.rm=TRUE)
                  deskrR[,"coefficient"] <- car::recode(deskrR[,"variable"], "'se.mean'='se';else='est'")
                  deskrR[,"parameter"]   <- gsub("se.mean","mean",deskrR[,"variable"])
                  deskrR   <- data.frame ( group = apply(deskrR[,allNam[["group"]],drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = allNam[["dependent"]], modus=modus, comparison = NA, deskrR[,c( "parameter", "coefficient", "value", allNam[["group"]])], stringsAsFactors = FALSE)
                  if(!is.null(allNam[["group.differences.by"]]))   {
                      if ( length (nCat ) >1 ) {
                           return(eatTools::facToChar(plyr::rbind.fill(deskrR,difs)))
                      }   else  {
                           return(eatTools::facToChar(deskrR))
                      }
                  }  else {
                      return(eatTools::facToChar(deskrR))
                  } }

jackknife.adjust.mean <- function (dat.i , allNam, na.rm, group.delimiter, type, repA, modus, scale , rscales, mse, rho, useEffectLiteR) {
          typeS<- car::recode(type, "'JK2'='JKn'")
          repl <- repA[ match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]]),]
          des  <- svrepdesign(data = dat.i[,c(allNam[["group"]], allNam[["dependent"]], allNam[["adjust"]])], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = 1, rscales = 1, repweights = repl[,-1, drop = FALSE], combined.weights = TRUE, mse = TRUE, rho=rho)
          if ( useEffectLiteR ) {
               ret <- withReplicates(des, funAdjustEL, allNam=allNam)
          }  else  {
               ret <- withReplicates(des, funAdjust, allNam=allNam)
          }
          rs   <- data.frame ( group = rep(attr(ret, "names"),2) , depVar = allNam[["dependent"]], modus = paste(modus, "survey", sep="__"), comparison = NA, parameter = "mean", coefficient = rep(c("est", "se"), each = nrow(as.data.frame (ret))), value = reshape2::melt(as.data.frame ( ret), measure.vars = colnames(as.data.frame ( ret)))[,"value"], rbind(do.call("rbind",  by(data=dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( x ) { x[1,allNam[["group"]], drop=FALSE]}, simplify = FALSE)),do.call("rbind",  by(data=dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( x ) { x[1,allNam[["group"]], drop=FALSE]}, simplify = FALSE))), stringsAsFactors=FALSE)
          return(rs)}
          
funAdjust <- function(w, data, allNam){                                         
       data[,allNam[["wgt"]]] <- w
       frml<- as.formula(paste0(allNam[["dependent"]]," ~ ", paste(allNam[["adjust"]], collapse = " + ")))
       if ( all(data[,allNam[["wgt"]]] == 1) )  {                               
            means <- by ( data = data, INDICES = data[,allNam[["group"]]], FUN = function ( gr ) {
                     reg <- lm(frml, data = gr)
                     cof1<- coef(reg)                                           
                     res1<- cof1[["(Intercept)"]]+ sum(unlist(lapply(names(cof1)[-1], FUN = function ( v ) { mean(data[,v]) * cof1[[v]]})))
                     return(res1)})
       }  else  {
            means <- by ( data = data, INDICES = data[,allNam[["group"]]], FUN = function ( gr ) {
                     do  <- paste0("lm(frml, data = gr, weights = ",allNam[["wgt"]],")")
                     reg <- eval(parse(text=do))
                     cof1<- coef(reg)
                     res1<- cof1[["(Intercept)"]]+ sum(unlist(lapply(names(cof1)[-1], FUN = function ( v ) { Hmisc::wtd.mean(data[,v], weights = data[,allNam[["wgt"]]]) * cof1[[v]]})))
                     return(res1)})
       }
       ret   <- as.vector(means)                                                
       names(ret) <- names(means)
       return(ret)}


funAdjustEL <- function(w, data, allNam){
       data[,allNam[["wgt"]]] <- w
       if ( all(data[,allNam[["wgt"]]] == 1) )  {                               
            res <- EffectLiteR::effectLite(y = allNam[["dependent"]], x = allNam[["group"]], z = allNam[["adjust"]], data = data, fixed.cell = TRUE, fixed.z = FALSE, homoscedasticity = FALSE, method="sem")
       }  else  {
            frml<- as.formula(paste0("~", allNam[["wgt"]]))
            res <- suppressMessages(EffectLiteR::effectLite(y = allNam[["dependent"]], x = allNam[["group"]], z = allNam[["adjust"]], data = data, fixed.cell = TRUE, fixed.z = FALSE, homoscedasticity = FALSE, weights = frml, method="sem"))
       }
       ret <- res@results@adjmeans[,"Estimate"]
       names(ret) <- names(table(data[,allNam[["group"]]]))
       return(ret)}  
       
conv.adjust.mean <- function ( dat.i, allNam, na.rm, group.delimiter, modus, useEffectLiteR) {
       if(isTRUE(useEffectLiteR)) {                                             
           if ( all(dat.i[,allNam[["wgt"]]] == 1) ) {                           
                res <- EffectLiteR::effectLite(y = allNam[["dependent"]], x = allNam[["group"]], z = allNam[["adjust"]], data = dat.i, fixed.cell = TRUE, fixed.z = FALSE, homoscedasticity = FALSE, method="sem")
           }  else  {
                st  <- paste("suppressMessages(EffectLiteR::effectLite(y = allNam[[\"dependent\"]], x = allNam[[\"group\"]], z = allNam[[\"adjust\"]], data = dat.i, fixed.cell = TRUE, fixed.z = FALSE, homoscedasticity = FALSE, weights = ~ ",allNam[["wgt"]],", method=\"sem\"))", sep="")
                res <- eval(parse( text=st))
           }
           vals <- reshape2::melt(res@results@adjmeans, measure.vars = c("Estimate", "SE"))[,"value"]
       }  else  {                                                               
           means <- do.call("rbind", by ( data = dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( gr ) {
                    frml<- paste0(allNam[["dependent"]], " ~ ", paste(allNam[["adjust"]], collapse = " + "))
                    if ( all(dat.i[,allNam[["wgt"]]] == 1) ) {                  
                         reg <- lm(as.formula(frml), data = gr)                 
                         x_m <- sapply(dat.i[,allNam[["adjust"]]], mean)        
                    }  else  {                                                  
                         reg <- eval(parse(text = paste0("lm(as.formula(frml), data = gr, weights = ",allNam[["wgt"]], ")")))
                         x_m <- unlist(lapply (allNam[["adjust"]], FUN = function (u) { Hmisc::wtd.mean(dat.i[,u], weights = dat.i[,allNam[["wgt"]]])}))
                    }                                                           
                    cof1<- coef(reg)
                    adj <- cof1[1] + sum(cof1[-1] * x_m)
                    pars<- c(cof1, x_m)
                    mat <- matrix(0, length(cof1) + length(x_m), length(cof1) + length(x_m))
                    mat[1:length(cof1), 1:length(cof1)] <- vcov(reg)
                    if ( all(dat.i[,allNam[["wgt"]]] == 1) ) {                  
                         mat[(length(cof1)+1):nrow(mat), (length(cof1)+1):ncol(mat)] <- var(gr[,allNam[["adjust"]]]) / (nrow(gr)-1)
                    }  else  {                                                  
                         mat[(length(cof1)+1):nrow(mat), (length(cof1)+1):ncol(mat)] <- cov.wt(gr[,allNam[["adjust"]], drop=FALSE], wt = gr[,allNam[["wgt"]]], cor = FALSE, center = TRUE)[["cov"]] / (nrow(gr)-1)
                    }
                    frm2<- paste0("~x1 + ", paste(paste("x", 2:length(cof1), sep=""), paste("x", (length(cof1)+1):(length(cof1)+length(x_m)), sep=""), collapse=" + ", sep="*"))
                    se  <- msm::deltamethod(as.formula(frm2), pars, mat)
                    return(data.frame ( mw = adj, se = se, stringsAsFactors = FALSE))}))
           vals  <- reshape2::melt(means, measure.vars = c("mw", "se"))[,"value"]
       }
       rs   <- data.frame ( group = rep(names(table(dat.i[,allNam[["group"]]])) , 2), depVar = allNam[["dependent"]], modus = modus, comparison = NA, parameter = "mean", coefficient = rep(c("est", "se"), each = length(vals)/2), value = vals, rbind(do.call("rbind",  by(data=dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( x ) { x[1,allNam[["group"]], drop=FALSE]}, simplify = FALSE)),do.call("rbind",  by(data=dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( x ) { x[1,allNam[["group"]], drop=FALSE]}, simplify = FALSE))), stringsAsFactors=FALSE)
       return(rs)}

jackknife.mean <- function (dat.i , allNam, na.rm, group.delimiter, type, repA, modus, scale, rscales, mse, rho) {
          typeS<- car::recode(type, "'JK2'='JKn'")
          repl <- repA[ match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]]),]
          des  <- svrepdesign(data = dat.i[,c(allNam[["group"]], allNam[["dependent"]])], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repl[,-1, drop = FALSE], combined.weights = TRUE, rho=rho)
          rets <- data.frame ( target = c("Ncases", "NcasesValid", "mean", "var"), FunctionToCall = c(NA,NA,"svymean","svyvar"), formelToCall = c("paste(\"~ \", \"N_weighted\",sep=\"\")","paste(\"~ \", \"N_weightedValid\",sep=\"\")","paste(\"~ \",allNam[[\"dependent\"]], sep = \"\")","paste(\"~ \",allNam[[\"dependent\"]], sep = \"\")"), naAction = c("FALSE","TRUE","na.rm","na.rm"), stringsAsFactors = FALSE)
          ret  <- apply(rets, 1, FUN = function ( toCall ) {                    
                  if (is.na(toCall[["FunctionToCall"]])) {                      
                      resL <- do.call("rbind", by(data=dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function (y){
                              if (toCall[["target"]] == "Ncases") {weg <- 0} else {weg <- length(which(is.na(y[,allNam[["dependent"]]])))}
                              r1 <- data.frame ( y[1,allNam[["group"]], drop=FALSE], parameter =toCall[["target"]], coefficient="est", value=nrow(y)-weg, stringsAsFactors = FALSE)
                              return(r1) }))
                  }  else  {
                      do   <- paste("svyby(formula = as.formula(",toCall[["formelToCall"]],"), by = as.formula(paste(\"~\", paste(allNam[[\"group\"]], collapse = \" + \"))), design = des, FUN = ",toCall[["FunctionToCall"]],",na.rm=",toCall[["naAction"]],", deff = FALSE, return.replicates = TRUE)",sep="")
                      res  <- suppressWarnings(eval(parse(text=do)))            
                      resL <- reshape2::melt( data = res, id.vars = allNam[["group"]], variable.name = "coefficient" , na.rm=TRUE)
                      stopifnot(length(table(resL[,"coefficient"])) == 2)
                      resL[,"coefficient"] <- car::recode(resL[,"coefficient"], "'se'='se'; else ='est'")
                      resL[,"parameter"]   <- toCall[["target"]]
                      attr(resL, "original") <- res
                  }
                  return(resL)})
          sds  <- do.call("rbind", by(data = dat.i, INDICES =  dat.i[,allNam[["group"]]], FUN = function (uu) {
                  namen   <- uu[1, allNam[["group"]], drop=FALSE]               
                  sub.rep <- repl[ match(uu[,allNam[["ID"]]], repl[,allNam[["ID"]]] ) ,  ]
                  des.uu  <- svrepdesign(data = uu[,c(allNam[["group"]], allNam[["dependent"]])], weights = uu[,allNam[["wgt"]]], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = sub.rep[,-1, drop = FALSE], combined.weights = TRUE, rho=rho)
                  var.uu  <- suppressWarnings(svyvar(x = as.formula(paste("~",allNam[["dependent"]],sep="")), design = des.uu, deff = FALSE, return.replicates = TRUE, na.rm = na.rm))
                  ret     <- data.frame(namen, est = as.numeric(sqrt(coef(var.uu))), se =  as.numeric(sqrt(vcov(var.uu)/(4*coef(var.uu)))), stringsAsFactors = FALSE )
                  return(ret)}) )
          sds  <- data.frame ( reshape2::melt(data = sds, id.vars = allNam[["group"]], variable.name = "coefficient" , na.rm=TRUE), parameter = "sd", stringsAsFactors = FALSE)
          resAl<- rbind(do.call("rbind",ret), sds)
          resAl<- data.frame ( group = apply(resAl[,allNam[["group"]],drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = allNam[["dependent"]], modus=paste(modus,"survey", sep="__"), comparison = NA, resAl[,c("parameter","coefficient","value",allNam[["group"]])] , stringsAsFactors = FALSE)
          if(!is.null(allNam[["group.differences.by"]]))   {
             nCat <- table(as.character(dat.i[,allNam[["group.differences.by"]]]))
             if ( length(nCat) < 2 ) {
                  warning("Grouping variable '", allNam[["group.differences.by"]], "' only has one category within imputation and/or nest. Group differences cannot be computed. Skip computation.")
             }  else  {
                m1   <- attr(ret[[ which(rets[,"target"] == "mean") ]], "original")
                sd   <- sds[which(sds[,"coefficient"] == "est"),]
                colnames(sd) <- car::recode(colnames(sd), "'value'='standardabweichung'")
                m    <- merge(m1, sd, by = allNam[["group"]], all = TRUE)
                stopifnot(nrow(m) == nrow(m1))
                m$comb.group <- apply(m, 1, FUN = function (ii) {eatTools::crop(paste( ii[allNam[["group"]]], collapse = "."))})
                repl1<- data.frame(t(attr(attr(ret[[which(rets[,"target"] == "mean")]], "original"), "replicates")), stringsAsFactors = FALSE )
                colnames(repl1) <- replCols <- paste("replNum", 1:ncol(repl1), sep="")
                repl1[,"comb.group"] <- rownames(repl1)
                m    <- merge(m, repl1, by = "comb.group" )
                m$all.group    <- 1
                res.group      <- tempR  <- setdiff(allNam[["group"]], allNam[["group.differences.by"]])
                if(length(res.group) == 0 ) {res.group <- "all.group"}
                kontraste      <- combinat::combn ( x = sort(unique(as.character(m[,allNam[["group.differences.by"]]]))), m = 2, simplify = FALSE)
                difs           <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                  ret <- do.call("rbind", lapply(kontraste, FUN = function ( k ) {
                                         if ( sum ( k %in% iii[,allNam[["group.differences.by"]]]) != length(k) ) {
                                              warning("Cannot compute contrasts for 'group.differences.by = ",allNam[["group.differences.by"]],"'.")
                                              return(NULL)                      
                                         } else {                               
                                              vgl.iii   <- iii[iii[,allNam[["group.differences.by"]]] %in% k ,]
                                              stopifnot ( nrow(vgl.iii) == 2 )
                                              true.diff <- diff(vgl.iii[,which(colnames(vgl.iii) %in% allNam[["dependent"]])])
                                              other.diffs <- apply(vgl.iii[,replCols], 2, diff)
                                              scumm     <- sapply(vgl.iii[,res.group,drop = FALSE], as.character)
                                              group     <- paste( paste( colnames(scumm), scumm[1,], sep="="), sep="", collapse = ", ")
                                              dummy     <- do.call("cbind", lapply ( allNam[["group"]], FUN = function ( gg ) {
                                                           ret <- data.frame ( paste ( unique(vgl.iii[,gg]), collapse = ".vs."))
                                                           colnames(ret) <- gg
                                                           return(ret)}))
                                              dif.iii   <- data.frame(dummy, group = group, vgl = paste(k, collapse = ".vs."), dif = true.diff, se =  sqrt(sum((true.diff - other.diffs)^2)), stringsAsFactors = FALSE )
                                              dif.iii   <- data.frame(dif.iii, es = dif.iii[["dif"]] / sqrt(0.5*sum(vgl.iii[,"standardabweichung"]^2)))
                                              return(dif.iii)
                                         } }))
                                  return(ret)}))
                difsL<- data.frame ( depVar = allNam[["dependent"]], reshape2::melt(data = difs, measure.vars = c("dif", "se", "es") , variable.name = "coefficient" , na.rm=TRUE), modus=paste(modus,"survey", sep="__"), parameter = "mean", stringsAsFactors = FALSE)
                difsL[,"coefficient"] <- car::recode(difsL[,"coefficient"], "'se'='se'; 'es'='es'; else = 'est'")
                difsL[,"comparison"]  <- "groupDiff"
                difsL[,"depVar"]      <- allNam[["dependent"]]
                difsL[,"group"] <- apply(difsL[,c("group","vgl")],1,FUN = function (z) {paste(z,collapse="____")})
                resAl<- rbind(resAl,difsL[,-match("vgl", colnames(difsL))])
             }
          }
          return(eatTools::facToChar(resAl)) }



buildString <- function(dat, allNam, refGrp, reihenfolge ) {
      strng <- paste0(rep(as.vector(unlist(by(data=dat, INDICES = dat[,allNam[["group"]]], FUN = function ( x ) {
               def <- c(x[1,allNam[["group"]]],refGrp[["groupValue"]])
               ref <- reihenfolge
               reih<- unlist(lapply(ref, FUN = function ( z ) {which(def %in% dat[,z])}))
               ret <- paste(def[reih], collapse="_")
               return(ret)}))), 2),giveRefgroup(refGrp))
      return(strng)}


giveRefgroup <- function ( refGrp) {
          if ( is.null(refGrp)) {return(".vs.wholeGroup")}
          ret <- paste0(".vs.", paste(apply(refGrp, MARGIN = 1, FUN = function ( zeile ) { zeile[["groupValue"]]}), collapse="_"))
          return(ret)}

jackknife.cov <- function (dat.i , allNam, na.rm, group.delimiter, type, repA, refGrp, scale , rscales , mse, rho, reihenfolge){
          typeS<- car::recode(type, "'JK2'='JKn'")
          repl <- repA[ match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]]),]
          des  <- svrepdesign(data = dat.i[,c(allNam[["group"]], allNam[["dependent"]])], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repl[,-1, drop = FALSE], combined.weights = TRUE, rho=rho)
          ret  <- withReplicates(des, groupVersusTotalMean, allNam=allNam)
          rs   <- data.frame ( group =  buildString(dat= dat.i,allNam=allNam, refGrp=refGrp, reihenfolge) , depVar = allNam[["dependent"]], modus = NA, comparison = "crossDiff", parameter = "mean", coefficient = rep(c("est", "se"), each = nrow(ret)), value = reshape2::melt(as.data.frame ( ret), measure.vars = colnames(as.data.frame ( ret)))[,"value"], rbind(do.call("rbind",  by(data=dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( x ) { x[1,allNam[["group"]], drop=FALSE]}, simplify = FALSE)),do.call("rbind",  by(data=dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( x ) { x[1,allNam[["group"]], drop=FALSE]}, simplify = FALSE))), stringsAsFactors=FALSE)
          return(rs)}


groupVersusTotalMean <- function(w, data, allNam){                              
       data[,allNam[["wgt"]]] <- w
       if ( all(data[,allNam[["wgt"]]] == 1) )  {
           ret <- by(data=data, INDICES = data[,allNam[["group"]]], FUN = function ( y ) { mean(y[,allNam[["dependent"]]])})
           gm  <- mean(data[,allNam[["dependent"]]])                            
       }  else  {
           ret <- by(data=data, INDICES = data[,allNam[["group"]]], FUN = function ( y ) { Hmisc::wtd.mean(y[,allNam[["dependent"]]], weights = y[,allNam[["wgt"]]])})
           gm  <- Hmisc::wtd.mean(data[,allNam[["dependent"]]], weights = data[,allNam[["wgt"]]])
       }                                                                        
       ret <- gm - ret
       return(ret)}

conv.cov <- function (dat.i, allNam, na.rm, group.delimiter, nBoot, refGrp, reihenfolge){
          covs<- boot::boot(data=dat.i, R = nBoot, statistic = function ( x, i) {groupVersusTotalMean(w = x[i,allNam[["wgt"]]], data = x[i,c(allNam[["group"]], allNam[["dependent"]])], allNam=allNam)})
          mns <- colMeans(covs$t)
          ses <- sapply(as.data.frame(covs$t), FUN = sd)                        
          rs  <- data.frame ( group =  buildString(dat= dat.i,allNam=allNam, refGrp=refGrp, reihenfolge) , depVar = allNam[["dependent"]], modus = NA, comparison = "crossDiff", parameter = "mean", coefficient = rep(c("est", "se"), each = length(mns)), value = c(mns, ses), rbind(do.call("rbind",  by(data=dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( x ) { x[1,allNam[["group"]], drop=FALSE]}, simplify = FALSE)),do.call("rbind",  by(data=dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( x ) { x[1,allNam[["group"]], drop=FALSE]}, simplify = FALSE))), stringsAsFactors=FALSE)
          return(rs)}
          
jackknife.glm <- function (dat.i , allNam, formula, forceSingularityTreatment, glmTransformation, na.rm, group.delimiter, type, repA, modus, useWec, scale, rscales, mse, rho, hetero, se_type, crossDiffSE.engine, stochasticGroupSizes, family) {
                 sub.ana <- by(data = dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function (sub.dat) {
                            nam    <- sub.dat[1,allNam[["group"]],drop=FALSE]
                            if ( allNam[["wgt"]] == "wgtOne") {
                                 glm.ii <- test <- glm(formula = formula, data = sub.dat, family = family)
                            }  else  {
                                 glm.ii <- test <- eval(parse(text = paste("glm(formula = formula, data = sub.dat, family = family, weights = ",allNam[["wgt"]],")",sep="")))
                            }
                            singular       <- names(glm.ii[["coefficients"]])[which(is.na(glm.ii[["coefficients"]]))]
                            if(!is.null(repA)) {
                                modus      <- paste(modus, "survey", sep="__")
                                typeS      <- car::recode(type, "'JK2'='JKn'")
                                design     <- svrepdesign(data = sub.dat[,c(allNam[["group"]], allNam[["independent"]], allNam[["dependent"]]) ], weights = sub.dat[,allNam[["wgt"]]], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(sub.dat[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                                if(length(singular) == 0 & isFALSE(forceSingularityTreatment) ) {
                                   glm.ii  <- suppressWarnings(svyglm(formula = formula, design = design, return.replicates = FALSE, family = family))
                                }
                            }
                            r.squared      <- data.frame ( r.squared = var(glm.ii$fitted.values)/var(glm.ii$y) , N = nrow(sub.dat) , N.valid = length(glm.ii$fitted.values) )
                            r.nagelkerke   <- fmsb::NagelkerkeR2(glm.ii)
                            summaryGlm     <- summary(glm.ii)
                            if( inherits(family, "family" )) {
                                if (  length(grep("binomial", eatTools::crop(capture.output(family)))) > 0 ) {
                                      res.bl <- data.frame ( group=paste(sub.dat[1,allNam[["group"]]], collapse=group.delimiter), depVar =allNam[["dependent"]],modus = modus, parameter = c(rep(c("Ncases","Nvalid",names(na.omit(glm.ii[["coefficients"]]))),2),"R2","R2nagel", "deviance", "null.deviance", "AIC", "df.residual", "df.null"),
                                                coefficient = c(rep(c("est","se"),each=2+length(names(na.omit(glm.ii[["coefficients"]])))),rep("est", 7)),
                                                value=c(r.squared[["N"]],r.squared[["N.valid"]],na.omit(glm.ii[["coefficients"]]),NA,NA,summaryGlm$coef[,2],r.squared[["r.squared"]],r.nagelkerke[["R2"]], test$deviance, test$null.deviance, test$aic, test$df.residual, test$df.null),sub.dat[1,allNam[["group"]], drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
                                }   else  {
                                      res.bl <- data.frame ( group=paste(sub.dat[1,allNam[["group"]]], collapse=group.delimiter), depVar =allNam[["dependent"]],modus = modus, parameter = c(rep(c("Ncases","Nvalid",names(na.omit(glm.ii[["coefficients"]]))),2),"R2","R2nagel"),
                                                coefficient = c(rep(c("est","se"),each=2+length(names(na.omit(glm.ii[["coefficients"]])))),rep("est", 2)),
                                                value=c(r.squared[["N"]],r.squared[["N.valid"]],na.omit(glm.ii[["coefficients"]]),NA,NA,summaryGlm$coef[,2],r.squared[["r.squared"]],r.nagelkerke[["R2"]]),sub.dat[1,allNam[["group"]], drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
                                }
                            }  else  {
                                res.bl <- data.frame ( group=paste(sub.dat[1,allNam[["group"]]], collapse=group.delimiter), depVar =allNam[["dependent"]],modus = modus, parameter = c(rep(c("Ncases","Nvalid",names(na.omit(glm.ii[["coefficients"]]))),2),"R2","R2nagel"),
                                          coefficient = c(rep(c("est","se"),each=2+length(names(na.omit(glm.ii[["coefficients"]])))),rep("est", 2)),
                                          value=c(r.squared[["N"]],r.squared[["N.valid"]],na.omit(glm.ii[["coefficients"]]),NA,NA,summaryGlm$coef[,2],r.squared[["r.squared"]],r.nagelkerke[["R2"]]),sub.dat[1,allNam[["group"]], drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
                            }
                            if(!is.null(repA) || isTRUE(useWec) ) {             
                                if(length(which(is.na(glm.ii[["coefficients"]]))) > 0 ) {
                                   message("Singularity problem in regression estimation for ", length(singular)," coefficient(s): ",paste(singular, collapse = ", "),". Try workaround ... ")
                                }
                                if(isTRUE(forceSingularityTreatment) && isFALSE(useWec)) {
                                   message("Compute coefficients in the expectation of singularities ... ")
                                }
                                if(length(singular) > 0 || forceSingularityTreatment == TRUE ) {
                                   stopifnot(length(as.character(formula)) == 3 )
                                   if ( isFALSE(useWec) ) {
                                       warning("Unidentified bug with Nagelkerkes r^2 in singularity treatment. No r^2 is computed.")
                                       if ( glmTransformation == "none" )  {resRoh <- withReplicates(design, getOutputIfSingular, allNam=allNam,  frml= formula, fam =family)}
                                       if ( glmTransformation == "sdY" )   {resRoh <- withReplicates(design, getOutputIfSingularT1, allNam=allNam,  frml= formula, fam =family)}
                                   }  else  {                                   
                                       if ( crossDiffSE.engine == "lavaan") {
                                            if(!is.null(repA) ) {               
                                                design1<- svrepdesign(data = dat.i[,c(allNam[["independent"]], allNam[["dependent"]]) ], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                                                ret    <- withReplicates(design1, funadjustLavaanWec, allNam=allNam)
                                                resRoh <- data.frame ( ret, stringsAsFactors=FALSE)
                                                rownames(resRoh) <- paste0(allNam[["independent"]], rownames(resRoh))
                                            }  else  {                          
                                                res    <- groupToTotalMeanComparisonLavaan ( d = dat.i, y.var = allNam[["dependent"]], group.var = allNam[["independent"]], weight.var=allNam[["wgt"]], heterogeneous=hetero, stchgrsz=stochasticGroupSizes, extended.results=FALSE, lavaan.syntax.path=NULL, run=TRUE, lavaan.summary.output=FALSE )
                                                resRoh <- data.frame ( theta = res[,"est"], SE = res[,"se"], stringsAsFactors=FALSE)
                                                rownames(resRoh) <- paste0(allNam[["independent"]], res[,"par"])
                                            }
                                       }  else  {                               
                                            contrasts(dat.i[,as.character(formula)[3]]) <- eatTools::contr.wec.weighted(dat.i[,as.character(formula)[3]], omitted=names(table(dat.i[,as.character(formula)[3]]))[1], weights = dat.i[,allNam[["wgt"]]])
                                            if(!is.null(repA) ) {
                                                design1<- svrepdesign(data = dat.i[,c(allNam[["independent"]], allNam[["dependent"]]) ], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                                                resRoh1<- data.frame ( withReplicates(design1, getOutputIfSingularWec, allNam=allNam, frml = formula), stringsAsFactors=FALSE)
                                                contrasts(dat.i[,as.character(formula)[3]]) <- eatTools::contr.wec.weighted(dat.i[,as.character(formula)[3]], omitted=names(table(dat.i[,as.character(formula)[3]]))[length(names(table(dat.i[,as.character(formula)[3]])))], weights =  dat.i[,allNam[["wgt"]]])
                                                design2<- svrepdesign(data = dat.i[,c(allNam[["independent"]], allNam[["dependent"]]) ], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                                                resRoh2<- data.frame ( withReplicates(design2, getOutputIfSingularWec, allNam=allNam, frml = formula), stringsAsFactors=FALSE)
                                                resRoh <- rbind(resRoh1, resRoh2)[which(!duplicated(c(row.names(resRoh1), row.names(resRoh2)))),]
                                            }  else  {                          
                                                res1   <- eval(parse(text=createCall ( hetero=hetero, allNam=allNam, formula=formula) ), envir=NULL)
                                                contrasts(dat.i[,as.character(formula)[3]]) <- eatTools::contr.wec.weighted(dat.i[,as.character(formula)[3]], omitted=names(table(dat.i[,as.character(formula)[3]]))[length(names(table(dat.i[,as.character(formula)[3]])))], weights =  dat.i[,allNam[["wgt"]]])
                                                res2   <- eval(parse(text=createCall ( hetero=hetero, allNam=allNam, formula=formula) ))
                                                rowNams<- c(rownames(summary(res1)[["coefficients"]]), rownames(summary(res2)[["coefficients"]]))
                                                weg    <- which(duplicated(rowNams))
                                                rowNams<- rowNams[-weg]         
                                                resRoh <- data.frame ( rbind ( summary(res1)[["coefficients"]], summary(res2)[["coefficients"]]), stringsAsFactors = FALSE)[-weg,1:2]
                                                colnames(resRoh) <- c("theta", "SE")
                                                resRoh <- rbind(resRoh, data.frame ( theta = c(summary(res1)[["r.squared"]], length(res1[["fitted.values"]])), SE = c(NA, NA), stringsAsFactors = FALSE))
                                                rownames(resRoh) <- c(rowNams, "R2", "Nvalid")
                                            }
                                        }
                                   }                                            
                                   index      <- which(nchar(rownames(resRoh)) == 0)
                                   if(length(index)>0) { for ( j in 1:length(index)) { rownames(resRoh)[index[j]] <- paste("dummyPar",j,sep="")}}
                                   res.bl     <- data.frame ( group=paste(sub.dat[1,allNam[["group"]]], collapse=group.delimiter), depVar =allNam[["dependent"]],modus = modus, parameter = rep(rownames(resRoh), 2), coefficient = c(rep("est", nrow(resRoh)), rep("se", nrow(resRoh))),
                                                 value = c(resRoh[,"theta"], resRoh[,"SE"]), sub.dat[1,allNam[["group"]], drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
                                   weg        <- intersect ( which(res.bl[,"coefficient"] == "se") , which(res.bl[,"value"] == 0) )
                                   if(length(weg)>0) { res.bl[weg,"value"] <- NA}
                                }
                            }
                            return(list ( res.bl=res.bl, glm.ii=glm.ii, nam=nam)) })
                 ori.glm <- lapply(sub.ana, FUN = function ( x ) { x[["glm.ii"]]})
                 nams    <- lapply(sub.ana, FUN = function ( x ) { x[["nam"]]})
                 sub.ana <- do.call("rbind", lapply(sub.ana, FUN = function ( x ) { x[["res.bl"]]}))
                 sub.ana[,"comparison"] <- NA
                 return(list(sub.ana=sub.ana, ori.glm=ori.glm, nams=nams)) }
                 
funadjustLavaanWec <- function(w, data, allNam){
       data[,allNam[["wgt"]]] <- w
       res  <- groupToTotalMeanComparisonLavaan ( d = data, y.var = allNam[["dependent"]], group.var = allNam[["independent"]], weight.var=allNam[["wgt"]], heterogeneous=TRUE, stchgrsz=FALSE, extended.results=FALSE, lavaan.syntax.path=NULL, run=TRUE, lavaan.summary.output=FALSE )
       ret  <- res[,"est"]
       names(ret) <- res[,"par"]
       return(ret)}

superSplitter <- function ( group=NULL, group.splits = length(group), group.differences.by = NULL, group.delimiter = "_" , dependent )  {
             group        <- as.list(group)
             names(group) <- unlist(group)
             if(max(group.splits)> length(group)) {group.splits[which(group.splits>length(group))] <- length(group)}
             group.splits <- unique(group.splits)
             superSplitti <- unlist(lapply(group.splits, FUN = function ( x ) {
                             spl <- combinat::combn(names(group),x)
                             if(inherits(spl, "matrix")) { spl <- as.list(data.frame(spl))} else {spl <- list(spl)}
                             spl <- unlist(lapply(spl, FUN = function ( y ) { paste(as.character(unlist(y)), collapse="________")}))
                             return(spl)}))
             superSplitti <- strsplit(superSplitti, "________")
             namen        <- unlist(lapply(superSplitti, FUN = function ( y ) { paste(y, collapse=group.delimiter)}))
             superSplitti <- lapply(superSplitti, FUN = function ( y ) {
                             if(!is.null(group.differences.by)) {if( group.differences.by %in% y ) { attr(y,"group.differences.by") <- group.differences.by}}
                             return(y)})
             names(superSplitti) <- namen
             return(superSplitti)}


dG <- function ( jk2.out , analyses = NULL, digits = 3, printDeviance, add ) {
            trend <- lapply(names(jk2.out[["resT"]]), FUN = function (tr) {
                     cat(paste("       Trend group: '", tr, "'.\n",sep=""))
                     splitData <- by ( data = jk2.out[["resT"]][[tr]], INDICES = jk2.out[["resT"]][[tr]][,c("group", "depVar")], FUN = function ( spl ) {return(spl)})
			               if(is.null(analyses)) {analyses <- 1:length(splitData)}
                     for ( i in analyses) {
                           spl    <- splitData[[i]]
                           weg1   <- which ( spl[,"parameter"] %in% c("Ncases","Nvalid","R2","R2nagel", "deviance", "df.null", "df.residual", "null.deviance", "AIC"))
                           weg2   <- grep("wholePopDiff",spl[,"parameter"])
                           weg3   <- grep("^p$", spl[,"coefficient"])
                           ret    <- reshape2::dcast(spl[-unique(c(weg1, weg2, weg3)),], parameter~coefficient)
                           ret[,"t.value"] <- ret[,"est"] / ret[,"se"]
                           df     <- spl[ spl[,"parameter"] == "Nvalid" & spl[,"coefficient"] == "est"  ,"value"] - nrow(ret)
                           ret[,"p.value"] <- 2*(1-pt( q = abs(ret[,"t.value"]), df = df ))
                           ret[,"sig"]     <- eatTools::num.to.cat(x = ret[,"p.value"], cut.points = c(0.001, 0.01, 0.05, 0.1), cat.values = c("***", "**", "*", ".", ""))
                           retNR  <- ret
                           ret    <- eatTools::roundDF(ret, digits=digits)
                           groupNamen <- setdiff(colnames(spl), c("group","depVar","modus", "parameter", "coefficient","value", "comparison"))
                           if ( length(groupNamen)>0) {
                                cat ( paste( "            groups: ", paste( groupNamen, unlist(lapply(spl[1,groupNamen], as.character)), sep=" = ", collapse = "; "),"\n",sep=""))
                           }
                           if ( length(add) >0 ) {
                                for ( jj in names(add)) {
                                     lz  <- 18 - nchar(jj)                      
                                     if ( lz < 0 ) {lz <- 0}
                                     cat(paste( paste(rep(" ", times = lz),collapse=""), jj, ": ", add[[jj]], "\n", sep=""))
                                }
                           }
                           cat ( paste( "dependent Variable: ", as.character(spl[1,"depVar"]), "\n \n", sep=""))
                           print(ret)
                           r2     <- spl[ spl[,"parameter"] == "R2" ,"value"]
                           r2nagel<- spl[ spl[,"parameter"] == "R2nagel" ,"value"]
                           cat(paste("\n            R-squared: ",round(r2[1],digits = digits),"; SE(R-squared): ",round(r2[2],digits = digits),"\n",sep=""))
                           cat(paste  ("Nagelkerkes R-squared: ",round(r2nagel[1],digits = digits),"; SE(Nagelkerkes R-squared): ",round(r2nagel[2],digits = digits),"\n",sep=""))
                           if ( isTRUE(printDeviance) ) {
                                for ( prms in c("deviance", "null.deviance", "df.null", "df.residual", "AIC")) {
                                      if ( prms %in% spl[,"parameter"] ) {
                                           cat(paste  ( paste(rep(" ", times = 21 - nchar(prms)),collapse=""),prms,": ", round ( spl[intersect(which(spl[,"parameter"] == prms), which(spl[,"coefficient"] == "est")),"value"],digits = digits), "\n", sep=""))
                                      }
                                }
                           }
                           nn     <- spl[ spl[,"parameter"] == "Nvalid" & spl[,"coefficient"] =="est" ,"value"]
                           cat(paste( round(nn, digits = 2), " observations and ",round(df,digits = 2), " degrees of freedom.",sep="")); cat("\n")
                           if(i != max(analyses)) { cat(paste0(paste(rep("-",times=66),collapse=""),"\n")) }
                     }
                     if ( tr != names(jk2.out[["resT"]])[length(names(jk2.out[["resT"]]))] ) {
                          cat(paste0(paste(rep("=",times=76),collapse=""),"\n"))
                     }
                     }) }


getOutputIfSingularWec <- function ( w, data, allNam, frml) {
                       data[,allNam[["wgt"]]] <- w
                       do  <- paste0("lm(formula = ",capture.output(frml)[1],", weights=",allNam[["wgt"]],", data=data)")
                       res <- eval(parse(text=do))
                       coefs <- c(res[["coefficients"]], R2 = summary(res)[["r.squared"]], Nvalid = length(res$fitted.values))
                       return(coefs)}

getOutputIfSingular <- function ( w, data, allNam, frml, fam ) {
                       data[,allNam[["wgt"]]] <- w
                       do  <- paste0("glm(formula = ",capture.output(frml),", weights=",allNam[["wgt"]],", data=data, family = fam)")
                       res <- eval(parse(text=do))
                       coefs <- na.omit(coef(res))
                       rnagel<- unlist(fmsb::NagelkerkeR2(res))
                       names(rnagel) <- c("Nvalid", "R2nagel")
                       rnagel<- rnagel[1]
                       coefs <- c(coefs, R2 = var(res$fitted.values)/var(res$y), rnagel)
                       return(coefs)}
                       
getOutputIfSingularT1<- function ( w, data, allNam, frml, fam ) {
                       data[,allNam[["wgt"]]] <- w
                       do  <- paste0("glm(formula = ",capture.output(frml),", weights=",allNam[["wgt"]],", data=data, family = fam)")
                       res <- eval(parse(text=do))
                       coefs <- na.omit(coef(res))
                       pred  <- sd ( res[["linear.predictors"]] ) +  (pi^2)/3
                       coefs <- coefs/pred
                       rnagel<- unlist(fmsb::NagelkerkeR2(res))
                       names(rnagel) <- c("Nvalid", "R2nagel")
                       rnagel<- rnagel[1]
                       coefs <- c(coefs, R2 = var(res$fitted.values)/var(res$y), rnagel)
                       return(coefs)}

clearTab <- function ( repTable.output, allNam , depVarOri, fc, toCall, datL) {
            if ( fc == "repTable" && toCall == "mean" ) {
                 stopifnot ( all(repTable.output[,"parameter"] %in% c("meanGroupDiff", "wholePopDiff") == FALSE))
                 jk2 <- repTable.output[which(repTable.output[,"parameter"] == "mean"),]
                 stopifnot(length(unique(jk2[,"depVar"]))==1)                   
                 Nc <- repTable.output[intersect(which(repTable.output[,"parameter"] == "Ncases"),  which(repTable.output[,"coefficient"] == "est")),]
                 if(!is.null(depVarOri)) {
                     prm <- datL[which(datL[,as.character(jk2[1,"depVar"])]==1),depVarOri]
                     stopifnot(length(unique(prm))==1)
                     jk2[,"depVar"]   <- depVarOri
                     Nc[,"depVar"]    <- depVarOri
                 }  else  {
                     prm <- "1"
                 }
                 jk2[,"parameter"] <- unique(prm)
                 jk2 <- rbind(jk2, Nc)
                 return(jk2)
            }  else  {
                 return(repTable.output)
            } }

chooseFunction <- function (datI, allNam, na.rm, group.delimiter,type, repA, modus, separate.missing.indicator ,
                  expected.values, probs,nBoot,bootMethod, formula, forceSingularityTreatment, glmTransformation, pb, toCall,
                  doJK, useWec, refGrp, scale, rscales, mse, rho, reihenfolge, hetero, se_type, useEffectLiteR, crossDiffSE.engine,
                  stochasticGroupSizes, correct, family) {
        pb$tick(); flush.console()
        if( toCall == "mean" ) {                                                
            if ( isTRUE(doJK) ) {
                 if (is.null(allNam[["adjust"]])) {
                     ana.i <- jackknife.mean (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, modus=modus, scale = scale, rscales = rscales, mse=mse, rho=rho)
                 }  else  { 
                     ana.i <- jackknife.adjust.mean (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, modus=modus, scale = scale, rscales = rscales, mse=mse, rho=rho, useEffectLiteR=useEffectLiteR)
                 }    
            }  else  {
                 if (is.null(allNam[["adjust"]])) {
                     ana.i <- conv.mean (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, modus=modus)
                 }  else  {
                     ana.i <- conv.adjust.mean (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, modus=modus, useEffectLiteR=useEffectLiteR)
                 }    
            }
        }
        if( toCall == "table" ) {
            if ( isTRUE(doJK) ) {
                 ana.i <- jackknife.table ( dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, separate.missing.indicator = separate.missing.indicator, expected.values=expected.values, modus=modus, scale = scale, rscales = rscales, mse=mse, rho=rho)
            }  else  {
                 ana.i <- conv.table (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, separate.missing.indicator = separate.missing.indicator, correct=correct, expected.values=expected.values, modus=modus)
            }
        }
        if( toCall == "quantile" ) {
            if ( isTRUE(doJK) ) {
                 ana.i <- jackknife.quantile (dat.i = datI, allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, probs=probs, modus=modus, scale = scale, rscales = rscales, mse=mse, rho=rho)
            }  else  {
                 ana.i <- conv.quantile (dat.i = datI, allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, probs=probs, nBoot=nBoot,bootMethod=bootMethod, modus=modus)
            }
        }
        if( toCall == "cov" ) {
            if ( isTRUE(doJK) ) {
                 ana.i <- jackknife.cov (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, refGrp=refGrp, scale = scale, rscales = rscales, mse=mse, rho=rho, reihenfolge=reihenfolge)
            }  else  {
                 ana.i <- conv.cov (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, nBoot=nBoot, refGrp=refGrp, reihenfolge=reihenfolge)
            }
        }
        glms <- grps <- nams <- NULL
        if( toCall == "glm" ) {
            doChek<- checkRegression ( dat = datI, allNam=allNam, useWec=useWec)
            ana.i <- jackknife.glm ( dat.i = datI , allNam=allNam, formula=formula, forceSingularityTreatment=forceSingularityTreatment, glmTransformation=glmTransformation, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, modus=modus, useWec=useWec, scale = scale, rscales = rscales, mse=mse, rho=rho, hetero=hetero, se_type=se_type, crossDiffSE.engine=crossDiffSE.engine, stochasticGroupSizes=stochasticGroupSizes, family=family)
            glms  <- ana.i[["ori.glm"]]
            grps  <- ana.i[["nams"]]
            nams  <- lapply(grps, FUN = function ( x ) { paste(as.character(unlist(x)), collapse=group.delimiter)})
            ana.i <- ana.i[["sub.ana"]]
        }
        ana.i <- data.frame ( ana.i, datI[1,c(allNam[["nest"]], allNam[["imp"]]),drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
        return(list(ana.i=ana.i, glms=glms, nams=nams, grps=grps))}


reconstructResultsStructureGlm <- function ( group, neu, grps, group.delimiter, pooled, allNam, modus, formula) {
        r2  <- lapply(neu[[group]], FUN = function ( x ) {var(x$fitted.values)/var(x$y)})
        r2n <- lapply(neu[[group]], FUN = function ( x ) {fmsb::NagelkerkeR2(x)[["R2"]]})
        Nval<- lapply(neu[[group]], FUN = function ( x ) {length(x$fitted.values)})
        mat <- which ( unlist(lapply(grps[[1]], FUN = function ( x ) { paste(as.character(unlist(x)), collapse = group.delimiter)})) == group)
        out <- summary(pooled[[group]])
        if ( "est" %in% colnames(out) ) {
             est <- "est"
        } else {
             stopifnot ( "estimate" %in% colnames(out))
             est <- "estimate"
        }
        if ( "se" %in% colnames(out) ) {
             se <- "se"
        } else {
             stopifnot ( "std.error" %in% colnames(out))
             se <- "std.error"
        }                                                                       
        prms<- all.vars(formula)[-1]
        col <- which(unlist(lapply(out, FUN = function ( co ) { length(unlist(lapply(prms, FUN = function ( l ) { grep(l, co)})))}))>0)
        stopifnot(length(col) <= 1)
        if ( length(col) == 1) { prm <- as.character(out[,col])}
        if ( length(col) == 0) { prm <- rownames(out)}
        ret <- data.frame ( group=group, depVar =allNam[["dependent"]],modus = modus, comparison = NA, parameter = c(rep(c("Ncases","Nvalid",prm),2),"R2","R2nagel"),
               coefficient = c(rep(c("est","se"),each=2+length(prm)),rep("est", 2)) , value= c(NA, min(unlist(Nval)),out[,est], NA, NA, out[,se], pool.R2(unlist(r2), unlist(Nval), quiet = TRUE )[["m.pooled"]],
               pool.R2(unlist(r2n), unlist(Nval), quiet = TRUE )[["m.pooled"]]), grps[[1]][[mat]], stringsAsFactors = FALSE, row.names = NULL)
        return(ret)}

checkData <- function ( sub.dat, allNam, toCall, separate.missing.indicator, na.rm) {
        if(!is.null(allNam[["PSU"]])) {
            nJkZones <- length(table(as.character(sub.dat[,allNam[["PSU"]]])))
            if(nJkZones<2)  {
               warning("Found group(s) with less than 2 PSUs. Please check your data!")
            }
        }                                                                       
        if( (toCall == "table" & isFALSE(separate.missing.indicator)) | (toCall %in% c("mean", "quantile", "glm") & isFALSE(na.rm) ) )  {
            nObserved <- length(which(is.na(sub.dat[, allNam[["dependent"]]])))
            if(nObserved>0) {                                                   
               if ( toCall %in% c("mean", "quantile", "glm") ) {
                    stop("Found unexpected missing data in dependent variable for at least one group. Please check your data or set 'na.rm==TRUE'.\n")
               }  else  {
                    warning("Found unexpected missing data in dependent variable for at least one group although 'separate.missing.indicator' was set to 'FALSE'.")
               }
            }
        }
        if ( toCall %in% c("mean", "quantile", "glm") & isFALSE(na.rm)) {       
            nMissing <- length(which(is.na(sub.dat[, allNam[["dependent"]]])))
            if(nMissing == nrow(sub.dat))  {
               stop("Some groups without any observed data. Please check your data!\n")
            }
        }
        return(sub.dat)}

checkNests <- function (x, allNam, toAppl, gr) {
        if(length(x[,allNam[["ID"]]]) != length(unique(x[,allNam[["ID"]]])))  {stop("ID variable '",allNam[["ID"]],"' is not unique within nests and imputations.")}
        if( length(toAppl[[gr]])>0) { ret <- lapply( toAppl[[gr]], FUN = function ( y ) {table(x[,y])}) } else {ret <- 1}
        return(list(ret=ret)) }


doSurveyAnalyses <- function (datL1, allNam, doJK, na.rm, group.delimiter, type, repA, modus, separate.missing.indicator, expected.values,
                    probs, nBoot,bootMethod, formula, forceSingularityTreatment, glmTransformation, toCall, poolMethod, useWec, refGrp, scale, rscales, mse, rho, reihenfolge, hetero, se_type, useEffectLiteR,
                    crossDiffSE.engine, stochasticGroupSizes, progress, correct, family, str1) {
        if(isTRUE(datL1[1,"isClear"])) {                                        
            nrep<- table(datL1[, c(allNam[["nest"]], allNam[["imp"]])])
            nrep<- prod(dim(nrep))
            cri1<- nrep > 4 & length(unique(datL1[,allNam[["ID"]]]))>2000 & doJK
            cri2<- FALSE
            if ( isFALSE(doJK)) {
                 cri2<- nrep > 9 & length(unique(datL1[,allNam[["ID"]]]))>5000
            }
            if ( progress && (isTRUE(cri1) | isTRUE(cri2)) ) {
                 pb  <- progress_bar$new( format = paste0(str1, " [:bar] :percent in :elapsed"), incomplete = " ", total = nrep, clear = FALSE, width= 75, show_after = 0.01)
            }  else  {
                 pb <- list()
                 pb$tick <- function (){return(NULL)}
            }
            ana <- do.call("rbind", by(data = datL1, INDICES = datL1[,allNam[["nest"]]], FUN = function ( datN ) {
                   anaI <- by(data = datN, INDICES = datN[,allNam[["imp"]]], FUN = chooseFunction, allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter,
                           type=type, repA=repA, modus=modus, separate.missing.indicator = separate.missing.indicator, expected.values=expected.values, probs=probs,
                           nBoot=nBoot,bootMethod=bootMethod, formula=formula, forceSingularityTreatment=forceSingularityTreatment, glmTransformation=glmTransformation, pb=pb,
                           toCall=toCall, doJK=doJK, useWec=useWec, refGrp=refGrp, scale = scale, rscales = rscales, mse=mse, rho=rho, reihenfolge=reihenfolge, hetero=hetero,
                           se_type=se_type, useEffectLiteR=useEffectLiteR, crossDiffSE.engine=crossDiffSE.engine, stochasticGroupSizes=stochasticGroupSizes,
                           correct=correct, family=family)
                   glms <- lapply(anaI, FUN = function ( x ) { x[["glms"]]})
                   nams <- lapply(anaI, FUN = function ( x ) { x[["nams"]]})
                   grps <- lapply(anaI, FUN = function ( x ) { x[["grps"]]})
                   if ( toCall == "glm" && poolMethod == "mice" && length(glms) > 1) {
                        aussen <- unlist(nams[[1]])
                        innen  <- names(nams)
                        for ( j in 1:length(glms)) { names(glms[[j]]) <- aussen }
                        neu    <- lapply(aussen, FUN = function ( a ) {lapply(innen, FUN = function ( i ) { glms[[i]][[a]]})})
                        pooled <- lapply(neu, FUN = function ( x ) { mice::pool(mice::as.mira(x)) } )
                        names(pooled) <- names(neu) <- aussen
                        anaI   <- do.call("rbind", lapply(aussen, FUN = reconstructResultsStructureGlm, neu=neu, grps=grps, group.delimiter=group.delimiter, pooled=pooled, allNam=allNam, modus=modus, formula=formula))
                   }  else  {                                                       
                        anaI <- do.call("rbind", lapply(anaI, FUN = function ( x ) { x[["ana.i"]]}))
                   }
                   return(anaI)}))
            if ( length(unique(ana[,"modus"])) >1 ) {
                 warning("Heterogeneous mode: '", paste(unique(ana[,"modus"]), collapse="', '"),"'")
            }
            mod <- unique(ana[,"modus"])
            if ( poolMethod == "mice" && toCall == "glm")  {
                 retList <- ana
            }  else  {
                 if( length(table(ana[,allNam[["imp"]]])) > 1 ) {
                     retList <- jk2.pool ( datLong = ana, allNam=allNam, forceSingularityTreatment = forceSingularityTreatment, modus=mod)
                 }  else  {
                     retList <- ana[,-match(c(allNam[["nest"]], allNam[["imp"]]), colnames(ana))]
                 }
            }
        }  else  {
            retList <- NULL
        }
        retList <- addSig (dat = retList, allNam = allNam)
        return(retList)}

checkEngine <- function ( engine , allNam, modus, toCall, type) {
          if (engine == "BIFIEsurvey") {
              if (!is.null(allNam[["nest"]]) ) {
                  message("Engine 'BIFIEsurvey' currently does not work for nested imputation. Set 'engine' to 'survey'.")
                  engine <- "survey"
              }
              if ( length(unlist(lapply(c("glm", "quantile"), FUN = function ( y ) {unlist(lapply(c(modus, toCall), FUN = function (w) { grep(y, w)}))}))) > 0 ) {
                  message("Engine 'BIFIEsurvey' currently does not work for regression models and quantiles. Set 'engine' to 'survey'.")
                  engine <- "survey"
              }
              if ( !type %in% c("JK2", "JK1", "NONE") ) {
                  message("Engine 'BIFIEsurvey' currently only works for jackknife 1 and jackknife 2. Set 'engine' to 'survey' due to type = '",type,"'.")
                  engine <- "survey"
              }
              if ( !is.null(allNam[["adjust"]])) {
                   message("Engine 'BIFIEsurvey' currently does not work for adjusted means. Set 'engine' to 'survey'.")
                   engine <- "survey"
              }                                                                 
              if ( toCall == "table" && !is.null(allNam[["group.differences.by"]]) ) {
                   message("Engine 'BIFIEsurvey' currently does not work for chi.square tests of group differences of frequency tables. Set 'engine' to 'survey'.")
                   engine <- "survey"
              }
          }
          return(engine)}

checkWecForUV <- function(dat, allNam) {
          if(!is.null(attr(dat, "wrapperForWec"))) {
             auchUV <- allNam[["independent"]]
          }  else  {
             auchUV <- NULL
          }
          return(auchUV)}

checkJK.arguments <- function(type, repWgt, PSU, repInd) {
          if ( type %in% c("JK1", "JK2")) {
               if ( is.null(repWgt)) {
                   if (is.null(PSU) ) {
                       stop("For type = '",type,"', 'PSU' must be defined if no replicate weights are specified (i.e., if 'repWgt' is NULL).")
                   }
                   if (type == "JK2" && is.null(repInd) ) {
                       stop("For type = '",type,"', 'repInd' must be defined if no replicate weights are specified (i.e., if 'repWgt' is NULL).")
                   }
               }
          }}

checkGroupVars <- function ( datL, allNam, auchUV) {
          if(!is.null(allNam[["group"]]) || !is.null(auchUV) ) {
             chk <- lapply(allNam[["group"]], FUN = function ( v ) { if ( !inherits(datL[,v],  c("factor", "character", "logical", "integer"))) {stop(paste0("Grouping variable '",v,"' must be of class 'factor', 'character', 'logical', or 'integer'.\n"))} })
             for ( gg in c(allNam[["group"]], auchUV) ) {                       ### levels der Gruppen duerfen keine "." oder "_" enthalten, falls cross differences berechnet werden sollen
                   if ( length(which(is.na(datL[,gg])))>0) { stop("Grouping variable '",gg,"' contains ",length(which(is.na(datL[,gg])))," missing values.")}
                   isImp<- TRUE                                                 
                   if ( is.null(allNam[["nest"]]) && is.null(allNam[["imp"]]) ) {
                        isImp <- FALSE                                          
                   }  else  {                                                   
                        imps  <- c(allNam[["nest"]], allNam[["imp"]])           
                        if ( length(imps) == 1) {
                             datL[,"g"] <- datL[,imps]
                        }  else  {
                             txt   <- paste0("datL |> tidyr::unite(\"g\", c(", paste(match(c(allNam[["nest"]],allNam[["imp"]]), colnames(datL)),collapse=", "),"), remove=FALSE)")
                             datL  <- eval(parse(text=txt))
                        }
                        isUniq<- eatGADS::checkUniqueness2(datL, varName=gg, idVar=allNam[["ID"]], impVar="g")
                        if(is.na(isUniq)) {
                            warning("'eatGADS::checkUniqueness2' returns NA. TRUE/FALSE is expected.")
                            isImp <- TRUE
                        }  else  {
                            if(isUniq) { isImp <- FALSE}
                        }
                   }
                   if ( isFALSE(isImp) ) {
                        if ( "g" %in% colnames(datL) ) {
                            d <- datL[which(datL[,"g"] == unique(datL[,"g"])[1]),]
                        }  else  {
                            d <- datL
                        }
                        chk2 <- lme4::isNested(d[,allNam[["ID"]]], d[,gg])
                   }  else  {
                        chk2 <- all(by(data = datL, INDICES = datL[,c(allNam[["nest"]], allNam[["imp"]])], FUN = function ( i ) { lme4::isNested(i[,allNam[["ID"]]], i[,gg])}))
                   }
                   if (isFALSE(chk2)) { warning("Grouping variable '",gg,"' is not nested within persons (variable '",allNam[["ID"]],"').") }
                   if ( inherits(datL[,gg], c("factor", "character")) && length(grep("\\.|_", datL[,gg])) > 0) {
                       message( "Levels of grouping variable '",gg, "' contain '.' and/or '_' which is not allowed. '.' and '_' will be deleted.")
                       if ( inherits ( datL[,gg], "factor")) {
                           levNew <- gsub("\\.|_", "", levels(datL[,gg]))
                           datL[,gg] <- factor(gsub("\\.|_", "", datL[,gg]), levels = levNew)
                       }  else  {
                           datL[,gg] <- gsub("\\.|_", "", datL[,gg])
                       }
                   }
             }
          }
          if(!is.null(allNam[["group"]]) | !is.null(allNam[["independent"]]) ) {
             for ( gg in c(allNam[["group"]], allNam[["independent"]]) ) {
                 if (inherits ( datL[,gg], "factor")) {                         
                     if ( any(table(datL[,gg]) == 0)) {
                          lev <- names(which(table(datL[,gg]) !=0))
                          nlv <- names(which(table(datL[,gg]) ==0))
                          message( "Delete level(s) '", paste(nlv, collapse="', '"), "' of grouping or independent variable '",gg,"' without any observations.")
                          datL[,gg] <- factor(as.character(datL[,gg]), levels =lev)
                     }
                 }
             }
          }
          if ( "g" %in% colnames(datL)) {datL <- datL[,-match("g", colnames(datL))]}
          return(datL)}

checkForAdjustmentAndLmer <- function(datL, allNam, groupWasNULL, formula.random, formula.fixed) {
          if(!is.null(allNam[["adjust"]]) || !is.null(formula.random) || !is.null(formula.fixed) ) {
             vars<- rbind(data.frame (type = rep("adjust", length(allNam[["adjust"]])),vars = allNam[["adjust"]] , stringsAsFactors = FALSE),
                          data.frame (type = rep("fixed", length(all.vars(formula.fixed))),vars = all.vars(formula.fixed) , stringsAsFactors = FALSE),
                          data.frame (type = rep("random", length(all.vars(formula.random))),vars = all.vars(formula.random) , stringsAsFactors = FALSE))
             for (v in unique(vars[,2])) {
                if ( inherits(datL[,v], c("logical")) ) {
                     message(paste0("Logical variable '",v,"' will be transformed into numeric."))
                     datL[,v] <- as.numeric(datL[,v])
                }
                if ( !inherits(datL[,v], c("numeric", "integer", "character", "factor"))) {stop(paste0("Adjusting variable '",v,"' must be of class 'numeric', 'integer', 'character', or 'factor'.\n")) }
                if ( length(which(is.na(datL[,v]))) >0 ) {stop(paste0("Adjusting variable '",v,"' has missing values."))}
             }
             vars[,"isFac"] <- sapply(vars[,2], FUN = function ( v ) {inherits(datL[,v], c("factor", "character"))})
             if (length(which(vars[,"isFac"]==TRUE))>0) {
                 for ( g in which(vars[,"isFac"]==TRUE)) {
                       if (vars[g,"type"] == "adjust") {
                            newFr <- model.matrix( as.formula (paste("~",vars[g,"vars"],sep="")), data = datL)[,-1,drop=FALSE]
                            message(paste("    Adjusting variable '",vars[g,"vars"],"' of class '",class(datL[,vars[g,"vars"]]),"' was converted to ",ncol(newFr)," indicator(s) with name(s) '",paste(colnames(newFr), collapse= "', '"), "'.",sep=""))
                            datL  <- data.frame(datL, newFr, stringsAsFactors=FALSE)
                            allNam[["adjust"]] <- c(setdiff(allNam[["adjust"]],vars[g,"vars"]), colnames(newFr))
                       }
                       if (vars[g,"type"] == "fixed") {
                            if ( !vars[g,"vars"] %in% extractFactorVarsFromFormula(formula.fixed) && vars[g,"vars"] %in% all.vars(formula.fixed) ) {
                                 stop(paste0("Variable '",vars[g,"vars"],"' of class '",class(datL[,vars[g,"vars"]]), "' should be specified with 'as.factor(",vars[g,"vars"],")' in the formula.fixed argument.\n  If '",vars[g,"vars"],"' should be modeled as numeric, please change variable class of '",vars[g,"vars"],"' in the data into numeric."))
                            }
                       }
                       if (vars[g,"type"] == "random") {
                            if ( !vars[g,"vars"] %in% extractFactorVarsFromFormula(formula.random) && vars[g,"vars"] %in% all.vars(formula.random) ) {
                                 warning(paste0("Variable '",vars[g,"vars"],"' of class '",class(datL[,vars[g,"vars"]]), "' should be specified with 'as.factor(",vars[g,"vars"],")' in the formula.random argument.\n  If '",vars[g,"vars"],"' should be modeled as numeric, please change variable class of '",vars[g,"vars"],"' in the data into numeric."))
                            }
                       }
                 }
             }
             if ( groupWasNULL && !is.null(allNam[["adjust"]])) {stop("When adjusted variables are defined, argument 'groups' must not be NULL.")}
             if ( length(allNam[["group"]])>1  && !is.null(allNam[["adjust"]]) ) {stop("When adjusted variables are defined, to date, only one grouping variable is allowed.")}
             if (!is.null(allNam[["group.differences.by"]]) && !is.null(allNam[["adjust"]])) {
                  message("Computation of group differences using 'group.differences.by' currently does not work for adjusted means.")
                  allNam[["group.differences.by"]] <- NULL
             }
          }
          return(list(allNam=allNam, datL=datL, formula.random=formula.random, formula.fixed=formula.fixed))}
          
checkNameConvention <- function( allNam) {
          na    <- c("isClear", "N_weightedValid", "N_weighted",  "wgtOne", "le")
          naGr  <- c("wholePop", "group", "depVar", "modus", "parameter", "coefficient", "value", "linkErr", "comparison", "sum", "trendvariable", "g", "le")
          naInd <- c("(Intercept)", "Ncases", "Nvalid", "R2",  "R2nagel", "linkErr")
          naGr1 <- which ( allNam[["group"]] %in% naGr )                        ### hier kuenftig besser: "verbotene" Variablennamen sollen automatisch umbenannt werden!
          if(length(naGr1)>0)  {stop(paste0("Following name(s) of grouping variables in data set are forbidden due to danger of confusion with result structure:\n     '", paste(allNam[["group"]][naGr1], collapse="', '"), "'\n  Please rename these variable(s) in the data set.\n"))}
          naInd1<- which ( allNam[["independent"]] %in% naInd )
          if(length(naInd1)>0)  {stop(paste0("Following name(s) of independent variables in data set are forbidden due to danger of confusion with result structure:\n     '", paste(allNam[["independent"]][naInd1], collapse="', '"), "'\n  Please rename these variable(s) in the data set.\n"))}
          na2   <- which ( unlist(allNam) %in% na )
          if(length(na2)>0)  {stop(paste0("Following variable name(s) in data set are forbidden due to danger of confusion with result structure:\n     '", paste(unlist(allNam)[na2], collapse="', '"), "'\n  Please rename these variable(s) in the data set.\n"))} }

setCrossDifferences <- function (cross.differences, allNam, group.splits) {
          if ( !is.list(cross.differences)) {
               if ( !is.logical (cross.differences)) {
                    message("Argument 'cross.differences' must be either logical or a list. Set 'cross.differences' to FALSE.")
                    cross.differences <- FALSE
               }
               if ( isTRUE(cross.differences)) {
                    if ( "wholeGroup" %in% allNam[["group"]] ) {
                          message("No groups defined. Set 'cross.differences' to FALSE.")
                          cross.differences <- FALSE
                    }  else  {
                          if ( length(group.splits)==1) {
                                message("Argument 'group.splits' was set to 1. No cross-level differences can be computed. Set 'cross.differences' to FALSE.")
                                cross.differences <- FALSE
                          }  else  {
                                cross.differences <- combinat::combn(x=group.splits,m=2, simplify = FALSE)
                          }
                    }
               }
          }  else  {                                                            
               if(!all(unlist(lapply(cross.differences, length)) == 2)) {stop("Each element in 'cross.differences' must be a vector of length 2.\n")}
               if(!all(unlist(lapply(cross.differences, inherits, what = c("integer","numeric"))))) {stop("Each element in 'cross.differences' must be a numerical vector.\n")}
               if(!all(unlist(lapply(cross.differences, FUN = function ( x ) { length(unique(x))})) ==2 ) ) {stop("Each element in 'cross.differences' must be a vector of 2 different numbers.\n")}
               vals <- unique(unlist(cross.differences))
               if(!all(vals %in% group.splits)) {stop("All numerical values in 'cross.differences' must be included in 'group.splits'.\n")}
               cross.differences <- lapply( cross.differences, sort )
               if ( length(cross.differences) != length(unique(cross.differences)) ) {
                    message("Some comparisons in 'cross.differences' are identical. Duplicated comparisons will be removed.")
                    cross.differences <- unique(cross.differences)
               }
          }
          allNam[["cross.differences"]] <- cross.differences
          return(allNam)}
          
createLoopStructure <- function(datL, allNam, verbose) {
          if( is.null(allNam[["imp"]]) )  { datL[,"imp"] <- 1; allNam[["imp"]] <- "imp" } else { stopifnot(length(allNam[["imp"]]) == 1 ); datL[,allNam[["imp"]]] <- as.character(datL[,allNam[["imp"]]])}
          if( is.null(allNam[["wgt"]]) )  { datL[,"wgtOne"] <- 1; allNam[["wgt"]] <- "wgtOne" } else {
              stopifnot(length(allNam[["wgt"]]) == 1 )
              if ( !inherits(datL[,allNam[["wgt"]]], c("numeric", "integer")) ) { stop ( paste("Error: 'wgt' variable '",allNam[["wgt"]],"' of class '",paste(class(datL[,allNam[["wgt"]]]),collapse = "', '"),"' has to be numeric.\n",sep="")) }
              isMis <- which(is.na(datL[,allNam[["wgt"]]]))
              isZero<- which ( datL[,allNam[["wgt"]]] == 0 )
              if(length(isMis)>0) { stop (paste ( "Error: Found ",length(isMis)," missing values in the weight variable '",allNam[["wgt"]],"'.\n",sep="")) }
              if(length(isZero)>0) { warning( "Found ",length(isZero)," zero weights in the weight variable '",allNam[["wgt"]],"'.") }
          }
          if(!is.null(allNam[["nest"]]))  {
              stopifnot(length(allNam[["nest"]]) == 1 )
              datL[,allNam[["nest"]]] <- as.character(datL[,allNam[["nest"]]])
              if(verbose){cat(paste("\nAssume nested structure with ", length(table(datL[,allNam[["nest"]]]))," nests and ",length(table(datL[,allNam[["imp"]]]))," imputations in each nest. This will result in ",length(table(datL[,allNam[["nest"]]]))," x ",length(table(datL[,allNam[["imp"]]]))," = ",length(table(datL[,allNam[["nest"]]]))*length(table(datL[,allNam[["imp"]]]))," imputation replicates.\n",sep=""))}
          }  else  { if(verbose){cat("\nAssume unnested structure with ",length(table(datL[,allNam[["imp"]]]))," imputations.\n",sep="")}}
          datL[,"isClear"] <- TRUE
          if( is.null(allNam[["nest"]]) ) { datL[,"nest"]  <- 1; allNam[["nest"]]  <- "nest" }
          if(!is.null(allNam[["group"]])) {                                     
              for ( jj in allNam[["group"]] )  { datL[,jj] <- as.character(datL[,jj]) }
          }
          return(list(datL = datL, allNam = allNam))}

assignReplicates <- function ( repWgt, allNam, datL, engine, type , progress, verbose) {
          if(!is.null(repWgt) ) {
              if ( !is.null(allNam[["PSU"]]) | !is.null(allNam[["repInd"]]) ) {
                    warning("Arguments 'PSU' and 'repInd' are expected to be NULL if replicate weights are already defined (via 'repWgt').\n    'PSU' and 'repInd' will be ignored.")
              }
          }
          if(!is.null(allNam[["repWgt"]]))  {
              repA <- data.frame ( datL[,allNam[["ID"]], drop=FALSE], datL[,allNam[["repWgt"]] ])
              repA <- repA[!duplicated(repA[,allNam[["ID"]]]),]
          }  else  {
              if(!is.null(allNam[["PSU"]]) && engine=="survey" )  {
                  repW <- datL[!duplicated(datL[,allNam[["ID"]]]),]
                  repA <- generate.replicates(dat = repW, ID = allNam[["ID"]], wgt = allNam[["wgt"]], PSU = allNam[["PSU"]], repInd = allNam[["repInd"]], type=type , progress=progress, verbose=verbose)
              }  else  { repA <- NULL}
          }
          return(repA)}

generate.replicates <- function ( dat, ID, wgt = NULL, PSU, repInd, type, progress, verbose )   {
          if(type %in% c("JK2", "BRR")) { stopifnot(length(PSU) == 1 & length(repInd) == 1 ) }
          if(type  == "JK1" ) { if(!is.null(repInd))  {
             cat("'repInd' is ignored for 'type = JK1'.\n")
             repInd <- NULL
          }  }
          allVars     <- list(ID = ID, wgt = wgt, PSU = PSU, repInd = repInd)
          all.Names   <- lapply(allVars, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = dat, variable=ii, warnIfMissing = FALSE)})
          dat.i       <- dat[,unlist(all.Names)]
          if(type %in% c("JK2", "BRR")) { if( !all( names(table(dat.i[,all.Names[["repInd"]]])) == c(0,1)) ) {stop("Only 0 and 1 are allowed for repInd variable.\n")} }
          zonen       <- names(table(as.character(dat.i[,all.Names[["PSU"]]]) ) )
          if ( verbose) { cat(paste("Create ",length(zonen)," replicate weights according to ",type," procedure.\n",sep=""))}
          if ( progress && nrow(dat)>2500 & length(zonen) > 50 ) {              
               pb     <- progress_bar$new( format = "         replicates [:bar] :percent in :elapsed", incomplete = " ", total = length(zonen), clear = FALSE, width= 75, show_after = 0.01)
          }
          missings    <- sapply(dat.i, FUN = function (ii) {length(which(is.na(ii)))})
          if(!all(missings == 0)) {
              mis.vars <- paste(names(missings)[which(missings != 0)], collapse = ", ")
              stop(paste("Found missing value(s) in variable(s) ", mis.vars,".\n",sep=""))
          }                                                                     
          reps <- data.frame ( lapply(zonen , FUN = function(ii) {              
                  if ( progress && nrow(dat)>2500 & length(zonen) > 50 ) { pb$tick() }
                  rep.ii <- dat.i[,all.Names[["wgt"]]]                          
                  if(type == "JK2")  { rep.ii[dat.i[,all.Names[["PSU"]]] == ii ] <- ifelse(dat.i[ dat.i[,all.Names[["PSU"]]] == ii ,all.Names[["repInd"]]] == 1, 0, 2 * rep.ii[dat.i[,all.Names[["PSU"]]] == ii ] ) }
                  if(type == "BRR")  { rep.ii <- ifelse(dat.i[ ,all.Names[["repInd"]]] == 1, 0, 2 * rep.ii ) }
                  if(type == "JK1")  {
                     rep.ii[ which ( dat.i[,all.Names[["PSU"]]] == ii) ] <- 0
                     rep.ii[ which ( dat.i[,all.Names[["PSU"]]] != ii) ] <- rep.ii[ which ( dat.i[,all.Names[["PSU"]]] != ii) ] *  ( sum(dat.i[,all.Names[["wgt"]]]) / sum (rep.ii))
                  }
                  return(rep.ii) }), stringsAsFactors = FALSE)
          colnames(reps) <- paste(all.Names[["wgt"]], 1:ncol(reps), sep="_")
          ret            <- data.frame(dat.i[,all.Names[["ID"]],drop=FALSE], reps, stringsAsFactors = FALSE)
          attr(ret, "n.replicates") <- length(zonen)
          return(ret) }

checkImpNest <- function (datL, doCheck, toAppl, gr, allNam, toCall, separate.missing.indicator, na.rm) {
          if(isTRUE(doCheck)) {                                     
             if ( length( toAppl[[gr]] ) > 1) {                     
                  crsTab <- table(datL[,toAppl[[gr]]])
                  if ( length(which(crsTab < 10 )) > 0 ) {
                       warning("Small number of observations in some combinations of grouping variables:\n   Recommend to remove these group(s).\n", eatTools::print_and_capture(crsTab, 3) )
                  }
             }
             impNes<- by(data = datL, INDICES = datL[, c(allNam[["nest"]], toAppl[[gr]]) ], FUN = function ( x ) { length(table(as.character(x[,allNam[["imp"]]])))}, simplify = FALSE)
             laenge<- which(sapply(impNes, length) == 0)
             if ( length(laenge ) > 0 ) {
                  warning(length(laenge), " combination(s) of groups without any observations. Analysis most probably will crash.")
             }
             impNes<- table(impNes[setdiff (1:length(impNes), laenge)])
             if(length(impNes) != 1 ) {warning("Number of imputations differ across nests and/or groups!\n", eatTools::print_and_capture(impNes, 3))}
             if(!is.null(allNam[["PSU"]]))  {
                  psuNes<- table ( by(data = datL, INDICES = datL[,allNam[["nest"]]], FUN = function ( x ) { length(table(as.character(x[,allNam[["PSU"]]])))}, simplify = FALSE) )
                  if(length(psuNes) != 1 ) {warning("Number of PSUs differ across nests!\n", eatTools::print_and_capture(psuNes, 3))}
             }
             impNes<- by(data = datL, INDICES = datL[, c(allNam[["nest"]], allNam[["imp"]]) ], FUN = checkNests, allNam=allNam, toAppl=toAppl, gr=gr, simplify = FALSE)
             impNes<- data.frame ( do.call("rbind", lapply(impNes, FUN = function ( x ) { unlist(lapply(x[["ret"]], FUN = length)) })) )
             if ( !all ( sapply(impNes, FUN = function ( x ) { length(table(x)) } ) == 1) ) { warning("Number of units in at least one group differs across imputations!")}
             datL  <- do.call("rbind", by(data = datL, INDICES = datL[,c( allNam[["group"]], allNam[["nest"]], allNam[["imp"]])], FUN = checkData, allNam=allNam, toCall=toCall, separate.missing.indicator=separate.missing.indicator, na.rm=na.rm))
             ok    <- table(datL[,"isClear"])
             if(length(ok) > 1 ) { message( ok[which(names(ok)=="FALSE")] , " of ", nrow(datL), " cases removed from analysis due to inconsistent data.") }
          }
          return(datL)}

prepExpecVal <- function (toCall, expected.values, separate.missing.indicator, allNam, datL) {
          if(toCall=="table") {
             misInd <- which(is.na(datL[,allNam[["dependent"]]]))
             if(isTRUE(separate.missing.indicator)) {
                if(length(misInd)>0) { datL[misInd,allNam[["dependent"]]] <- "<NA>"}
             }  else {
                if(length(misInd)>0) {
                   warning("No seperate missing categorie was chosen. ", length(misInd), " missings were found anyhow for ",allNam[["dependent"]],". Missings will be deleted from the data.")
                   if(length(misInd) == nrow(datL)) {stop()}
                   datL <- datL[-misInd,]
                }
             }
             expected.values <- sort(unique(c(expected.values, names(table(datL[,allNam[["dependent"]]])))))
          }
          return(list(datL=datL, expected.values=expected.values))}

createInfoString <- function ( ai, toCall, gr, toAppl) {
          nr <- ifelse(is.null(ai), 1, match(gr, names(toAppl)))
          st <- paste0(paste(rep(" ", times = 8-nchar(toCall)),collapse=""), toCall, " analysis ",nr)
          return(st)}

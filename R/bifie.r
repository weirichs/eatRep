doBifieAnalyses <- function (dat.i, allNam, na.rm, group.delimiter,separate.missing.indicator, expected.values, probs, formula, glmTransformation, toCall, modus, type, verbose, L1wgt, L2wgt,formula.fixed, formula.random){
      dat.i<- eatTools::facToChar(dat.i[,intersect(unlist(allNam), colnames(dat.i))], from = "character", to = "factor")
      dat.g<- eatGADS::import_DF(dat.i, checkVarNames = FALSE)                  
      dat2 <- eatGADS::extractData(dat.g, convertLabels = "numeric")            
      dat2 <- sortBifie(dat = dat2, toCall=toCall, allNam=allNam)
      labsD<- dat.g[["labels"]][which(dat.g[["labels"]][,"varName"] == allNam[["dependent"]]),]
      datL <- by ( data = dat2, INDICES = dat2[,allNam[["imp"]]], FUN = function ( imp.i ) { return(imp.i)})
      jkt  <- car::recode(type, "'JK2'='JK_TIMSS'; 'JK1'='JK_GROUP'")
      txt  <- capture.output(bo   <- BIFIE.data.jack( data= datL,  wgt = allNam[["wgt"]], jktype=jkt , jkzone = allNam[["PSU"]], jkrep = allNam[["repInd"]], cdata=FALSE))
      if ( isTRUE(verbose)) { cat("\n"); print(bo)}
      attributes(allNam[["group"]]) <- NULL                                     
      if ( toCall == "mean") {
           resML <- bifieMean(bifie.obj = bo, allNam=allNam, dat.g=dat.g, labsD=labsD, toCall=toCall, dat.i=dat.i, modus=modus)
      }
      if ( toCall == "table") {
           resML <- bifieTable(bifie.obj = bo, allNam=allNam, dat.g=dat.g, labsD=labsD, toCall=toCall, dat.i=dat.i, modus=modus)
      }
      if ( toCall == "lmer") {
           resML <- bifieLmer(bifie.obj=bo, allNam=allNam, dat.g=dat.g, labsD=labsD, modus=modus, formula.fixed=formula.fixed, formula.random=formula.random)
      }
      return(resML)}


sortBifie <- function(dat, toCall, allNam) {
      if ( toCall != "lmer") {
           dat <- dat[order(dat[,allNam[["ID"]]]), ]                            
      }  else  {
           dat <- dat[order(dat[,allNam[["clusters"]]]), ]                      
      }
      return(dat)}

aufbOut <- function (resM, allNam, dat.g, labsD) {
      if(length(allNam[["group"]])==1) {
           cols <- grep("groupva", colnames(resM[["stat"]]), value=TRUE)
           stopifnot(length(cols) ==2)
           altn <- data.frame ( alt = cols, neu = paste0(cols, "1"), stringsAsFactors = FALSE)
           recs <- paste("'",altn[,"alt"] , "' = '" , altn[,"neu"],"'",sep="", collapse="; ")
           colnames(resM[["stat"]]) <- car::recode(colnames(resM[["stat"]]), recs)
      }
      if(length(allNam[["group"]])>0) {
           for ( i in 1:length(allNam[["group"]])) {
                 labsG<- dat.g[["labels"]][which(dat.g[["labels"]][,"varName"] == allNam[["group"]][i]),]
                 recs <- paste("'",labsG[,"value"] , "' = '" , labsG[,"valLabel"],"'",sep="", collapse="; ")
                 resM[["stat"]][,paste0("groupval", i)] <- car::recode (resM[["stat"]][,paste0("groupval", i)], recs)
           }
      }
      if(!is.na(labsD[1,"value"])) {                                            
           recs <- paste("'",labsD[,"value"] , "' = '" , labsD[,"valLabel"],"'",sep="", collapse="; ")
           resM[["stat"]][,"varval"] <- car::recode(resM[["stat"]][,"varval"], recs)
      }                                                                         
      return(resM)}

aufbNoMultilevel <- function (resM, mv, toCall, allNam, dat.i, modus, grp) {
          idv  <- c(grep("varval", colnames(resM[["stat"]]), value=TRUE), grep("groupval", colnames(resM[["stat"]]), value=TRUE))
          resML<- eatTools::facToChar(reshape2::melt(data=resM[["stat"]], id.vars=idv, measure.vars=mv,  na.rm=TRUE))
          resML[setdiff(1:nrow(resML), grep("_", resML[,"variable"])),"variable"] <- paste(resML[setdiff(1:nrow(resML), grep("_", resML[,"variable"])),"variable"], "est", sep="_")
          resML<- data.frame ( resML, reshape2::colsplit(string = resML[,"variable"], pattern="_", names = c("parameter", "coefficient")),stringsAsFactors = FALSE)
          if ( toCall == "table") {
               resML <- resML[which(resML[,"parameter"] == "perc"),]
               resML[,"parameter"] <- resML[,"varval"]
               resML[,"variable"]  <- resML[,"varval"] <- NULL                  
          }  else {
               resML[,"parameter"] <- car::recode(resML[,"parameter"], "'M'='mean'; 'SD'='sd'; 'Nweight'='NcasesValid'")
               resML[,"variable"]  <- NULL
          }
          resML[,"coefficient"] <- car::recode(resML[,"coefficient"], "'SE'='se'")
          recs <- paste("'",grep("groupval", colnames(resML), value=TRUE) , "' = '" , allNam[["group"]],"'",sep="", collapse="; ")
          colnames(resML) <- car::recode(colnames(resML), recs)
          if ( toCall == "table") {                                             
               Ns   <- do.call("rbind", by(dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function (y) {data.frame ( y[1,allNam[["group"]], drop=FALSE], parameter = "Ncases", coefficient="est", value=length(unique(y[,allNam[["ID"]]])), stringsAsFactors = FALSE) }))
               resML<- plyr::rbind.fill(resML, Ns)
          }
          resML[,"modus"] <- paste(modus, "BIFIEsurvey", sep="__")
          resML[,"depVar"]<- allNam[["dependent"]]
          resML[,"comparison"] <- NA
          if ( length(allNam[["group"]]) > 1) {
               resML[,"group"] <- apply(resML[,allNam[["group"]]], MARGIN = 1, FUN = function ( z ) {paste(z, collapse="_")})
          }
          if ( length(allNam[["group"]]) == 1) { resML[,"group"] <- resML[,allNam[["group"]]]}
          resML<- rbind(resML, grp)
          return(resML)}

bifieMean <- function ( bifie.obj, allNam, dat.g, labsD, toCall, dat.i, modus) {
      txt  <- capture.output(resM <- BIFIE.univar( BIFIEobj=bifie.obj , vars = allNam[["dependent"]], group=allNam[["group"]] ))
      mv   <- c("Nweight", "Ncases", "M", "M_SE", "M_p", "SD", "SD_SE", "SD_p")
      grp  <- computeGroupDifferences(resM=resM, allNam=allNam, dat.g=dat.g, modus=modus)
      resM <- aufbOut(resM=resM, allNam=allNam, dat.g=dat.g, labsD=labsD)
      resML<- aufbNoMultilevel(resM=resM, mv=mv, toCall=toCall, allNam=allNam, dat.i=dat.i, modus=modus, grp=grp)
      return(resML)  }

bifieTable <- function ( bifie.obj, allNam, dat.g, labsD, toCall, dat.i, modus) {
      txt  <- capture.output(resM <- BIFIE.freq( BIFIEobj=bifie.obj , vars = allNam[["dependent"]], group=allNam[["group"]] ))
      resM[["stat"]][,"perc_p"] <- 2*(1-pnorm(resM[["stat"]][,"perc"] / resM[["stat"]][,"perc_SE"]))
      mv   <- c("Nweight", "Ncases", "perc", "perc_SE", "perc_p")               
      stopifnot(is.null(allNam[["group.differences.by"]]))                      
      resM <- aufbOut(resM=resM, allNam=allNam, dat.g=dat.g, labsD=labsD)
      resML<- aufbNoMultilevel(resM=resM, mv=mv, toCall=toCall, allNam=allNam, dat.i=dat.i, modus=modus, grp=NULL)
      return(resML)  }
      
bifieLmer <- function ( bifie.obj, allNam, dat.g, labsD, modus, formula.fixed, formula.random) {
      resM <- BIFIE.twolevelreg( BIFIEobj=bifie.obj, dep=allNam[["dependent"]], formula.fixed=formula.fixed, formula.random=formula.random, idcluster=allNam[["clusters"]], wgtlevel2=allNam[["L2wgt"]], wgtlevel1=allNam[["L1wgt"]], group = allNam[["group"]])
      resM <- aufbOut(resM=resM, allNam=allNam, dat.g=dat.g, labsD=labsD)       
      resML<- aufbMultilevel(resM=resM, allNam=allNam, modus=modus)
      return(resML)}

computeGroupDifferences <- function(resM, allNam, dat.g, modus){
      if(is.null(allNam[["group.differences.by"]])) {return(NULL)} else {       
           liste<- data.frame ( merge(resM[["stat_M"]],resM[["stat_SD"]][,c(grep("^groupva", colnames(resM[["stat_SD"]]), value=TRUE), "SD", "SD_SE")],  by = grep("^groupva", colnames(resM[["stat_SD"]]), value=TRUE), all=TRUE, sort=FALSE), dp = resM[["parnames"]], groupvar0 = "wholePop", groupval0 = 0, stringsAsFactors = FALSE)
           if ( length(allNam[["group.differences.by"]]) == length(allNam[["group"]])) {
               colnames(liste) <- car::recode(colnames(liste), "'groupvar'='groupvar1'; 'groupval'='groupval1'")
           }                                                                    
           col  <- colnames(liste)[which(sapply(eatTools::facToChar(liste[1,]), FUN = function ( x ) { x == allNam[["group.differences.by"]] }))]
           col  <- c(eatTools::removeNumeric(col), eatTools::removeNonNumeric(col))
           res  <- setdiff(grep("^groupval", colnames(liste), value=TRUE), paste0("groupval", col[2]))
           grp  <- do.call("rbind", by(data=liste, INDICES = liste[,res], FUN = function ( x ) {
                   comb <- data.frame ( combinat::combn(x=x[,"dp"], m=2), stringsAsFactors = FALSE)
                   diffs<- do.call("rbind", lapply(comb, FUN = function ( y ) { 
                           dp   <- eval(parse(text=paste("list ( \"groupDiff\" =~ 0 + I(",y[1],"-",y[2],"))")))
                           resMd<- BIFIE.derivedParameters( resM, derived.parameters=dp )
                           if ( length ( res ) == 1) {                          
                                rg  <- "all.group=1"                            
                           }  else  {                                           
                                rg  <- setdiff(eatTools::removeNonNumeric(res),0)
                                nam <- lapply(rg, FUN = function ( z ) {
                                       nam <- liste[1,grep(paste0("groupvar", z), colnames(liste))]
                                       val <- x[1,paste0("groupval", z)]
                                       ret <- data.frame ( nam = nam, wert = val, toRec = dat.g[["labels"]][intersect(which(dat.g[["labels"]][,"varName"] == nam), which(dat.g[["labels"]][,"value"] == val)),"valLabel"], stringsAsFactors = FALSE)
                                       ret <- paste0(ret[["nam"]],"=",ret[["toRec"]])
                                       return(ret)})
                                rg  <- paste(unlist(nam), collapse=", ")
                           }                                                    
                           vs   <- x[which(x[,"dp"] %in% y),paste(c("groupval",col[2]),collapse="")]
                           vs   <- paste(dat.g[["labels"]][intersect(which(dat.g[["labels"]][,"varName"] %in% allNam[["group.differences.by"]]), which(dat.g[["labels"]][,"value"] %in% vs)),"valLabel"], collapse=".vs.")
                           es   <- resMd[["stat"]][["coef"]] / sqrt(0.5*sum(x[which(x[,"dp"] %in% y),"SD"]^2))
                           add  <- x[1,gsub("val", "var", setdiff(res, "groupval0")), drop=FALSE]
                           if ( ncol (add)>0) {
                              val <- x[1,setdiff(res, "groupval0"),drop=FALSE]
                              colnames(val) <- unlist(add)
                              val <- data.frame(t(sapply(colnames(val), FUN = function (v) {dat.g[["labels"]][intersect(which(dat.g[["labels"]][,"varName"] == v), which(dat.g[["labels"]][,"value"] == val[[v]])),"valLabel"]})), stringsAsFactors = FALSE)
                              add <- data.frame ( val, data.frame ( v1=vs, stringsAsFactors = FALSE))
                              colnames(add)[ncol(add)] <- allNam[["group.differences.by"]]
                           } else {
                              add <- data.frame ( v1=vs, stringsAsFactors = FALSE)
                              colnames(add) <- allNam[["group.differences.by"]]
                           }                                                    
                           ret  <- data.frame ( group = paste(rg, vs, sep="___"), depVar = allNam[["dependent"]], modus = paste(modus,"BIFIEsurvey", sep="__"),  comparison = "groupDiff", parameter = "mean", coefficient = c("est", "se", "p", "es"), value = c( (-1) * resMd[["stat"]][["coef"]],resMd[["stat"]][["se"]],resMd[["stat"]][["p"]], (-1)*es),add, stringsAsFactors = FALSE)
                           return(ret)}))
                   return(diffs)}))
      }
      return(grp)}

aufbMultilevel <- function(resM, allNam, modus) {
      if(!is.null(allNam[["group"]])) {
          colnames(resM[["stat"]]) <- car::recode(colnames(resM[["stat"]]), paste0("'groupval1'='",allNam[["group"]],"'"))
          resM[["stat"]][,"group"] <- resM[["stat"]][,allNam[["group"]]]
      }  else  {
          resM[["stat"]][,"group"] <- "wholeGroup"
      }                                                                         
      resL <- eatTools::facToChar(reshape2::melt(data=resM[["stat"]], id.vars=c("parameter", "group", allNam[["group"]]), measure.vars=c("est", "SE", "p"),  na.rm=TRUE))
      if(!is.null(allNam[["group"]])) {
          resML<- data.frame ( group=unique(resM[["stat"]][,"group"]),depVar = allNam[["dependent"]], modus=modus, parameter="Nvalid", coefficient="est",value=resM[["Npers"]],gruppenname = unique(resM[["stat"]][,allNam[["group"]]]), stringsAsFactors = FALSE)
          rows <- grep("^beta_", resL[,"parameter"])
          resML<- rbind(resML,data.frame ( group=resL[rows,"group"],depVar = allNam[["dependent"]], modus=modus, parameter=eatTools::removePattern(resL[rows,"parameter"],"beta_"), coefficient=tolower(resL[rows,"variable"]),value=resL[rows,"value"], gruppenname = resL[rows,allNam[["group"]]], stringsAsFactors = FALSE))
          rows <- grep("^ICC_|^R2", resL[,"parameter"])
          resML<- rbind(resML,data.frame ( group=resL[rows,"group"],depVar = allNam[["dependent"]], modus=modus, parameter=resL[rows,"parameter"], coefficient=tolower(resL[rows,"variable"]),value=resL[rows,"value"], gruppenname = resL[rows,allNam[["group"]]], stringsAsFactors = FALSE))
          colnames(resML) <- car::recode(colnames(resML), paste0("'gruppenname'='",allNam[["group"]],"'"))
      }  else  {
          resML<- data.frame ( group="wholeGroup", depVar = allNam[["dependent"]], modus=modus, parameter="Nvalid",coefficient="est",value=resM[["N"]], stringsAsFactors = FALSE)
          rows <- grep("^beta_", resL[,"parameter"])
          resML<- rbind(resML,data.frame ( group=resL[rows,"group"],depVar = allNam[["dependent"]], modus=modus, parameter=eatTools::removePattern(resL[rows,"parameter"],"beta_"), coefficient=tolower(resL[rows,"variable"]),value=resL[rows,"value"],  stringsAsFactors = FALSE))
          rows <- grep("^ICC_|^R2", resL[,"parameter"])
          resML<- rbind(resML,data.frame ( group=resL[rows,"group"],depVar = allNam[["dependent"]], modus=modus, parameter=resL[rows,"parameter"], coefficient=tolower(resL[rows,"variable"]),value=resL[rows,"value"], stringsAsFactors = FALSE))
      }
      return(resML)}

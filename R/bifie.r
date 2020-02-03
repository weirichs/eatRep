doBifieAnalyses <- function (dat.i, allNam, na.rm, group.delimiter,separate.missing.indicator, expected.values, probs, formula, glmTransformation, toCall, modus){
      dat.i<- facToChar(dat.i, from = "character", to = "factor")
      dat.g<- import_DF(dat.i, checkVarNames = FALSE)                           ### dies hier geschieht alles, um die Variablen numerisch zu machen, BIFIEsurvey will es so
      dat2 <- extractData(dat.g, convertLabels = "numeric")
      dat2 <- dat2[order(dat2[,allNam[["ID"]]]), ]
      labsD<- dat.g[["labels"]][which(dat.g[["labels"]][,"varName"] == allNam[["dependent"]]),]
   ### in Listenformat bringen, to do: genestete Imputationen?
      datL <- by ( data = dat2, INDICES = dat2[,allNam[["imp"]]], FUN = function ( imp.i ) { return(imp.i)})
      txt  <- capture.output(bo   <- BIFIE.data.jack( data= datL,  wgt = allNam[["wgt"]], jktype="JK_TIMSS" , jkzone = allNam[["PSU"]], jkrep = allNam[["repInd"]], cdata=FALSE ))
      if ( toCall == "mean") {
           txt  <- capture.output(resM <- BIFIE.univar( BIFIEobj=bo , vars = allNam[["dependent"]], group=allNam[["group"]] ))
           mv   <- c("Nweight", "Ncases", "M", "M_SE", "M_p", "SD", "SD_SE", "SD_p")
      }                                                                         ### obere Zeile: dies sind die "measure.vars" fuers reshapen
      if ( toCall == "table") {
           txt  <- capture.output(resM <- BIFIE.freq( BIFIEobj=bo , vars = allNam[["dependent"]], group=allNam[["group"]] ))
           resM[["stat"]][,"perc_p"] <- 2*(1-pnorm(resM[["stat"]][,"perc"] / resM[["stat"]][,"perc_SE"]))
           mv   <- c("Nweight", "Ncases", "perc", "perc_SE", "perc_p")          ### measure.vars fuers reshapen
      }                                                                         ### jetzt Gruppendifferenzen ueber den BIFIE-hotfix
      if(!is.null(allNam[["group.differences.by"]])) {                          ### ueberoberunschoen! die 'parnames' haben nur fuer die BIFIE.univar-Funktion komische krude Namen, deswegen ist das hier alles komplizierter als es die Hilfeseite von BIFIE.by andeutet ... und man muss UNBEDINGT mit sort=FALSE mergen, sonst stimmt die Reihenfolge in 'resM' nicht mit der in 'redMd' ueberein!
           liste<- data.frame ( merge(resM[["stat_M"]],resM[["stat_SD"]][,c(grep("^groupva", colnames(resM[["stat_SD"]]), value=TRUE), "SD", "SD_SE")],  by = grep("^groupva", colnames(resM[["stat_SD"]]), value=TRUE), all=TRUE, sort=FALSE), dp = resM[["parnames"]], groupvar0 = "wholePop", groupval0 = 0, stringsAsFactors = FALSE)
           if ( length(allNam[["group.differences.by"]]) == length(allNam[["group"]])) {
               colnames(liste) <- recode(colnames(liste), "'groupvar'='groupvar1'; 'groupval'='groupval1'")
           }                                                                    ### boah, ich werd wahnsinnig ...
           col  <- colnames(liste)[which(sapply(facToChar(liste[1,]), FUN = function ( x ) { x == allNam[["group.differences.by"]] }))]
           col  <- c(removeNumeric(col), removeNonNumeric(col))
           res  <- setdiff(grep("^groupval", colnames(liste), value=TRUE), paste0("groupval", col[2]))
           grp  <- do.call("rbind", by(data=liste, INDICES = liste[,res], FUN = function ( x ) {
                   comb <- data.frame ( combn(x=x[,"dp"], m=2), stringsAsFactors = FALSE)
                   diffs<- do.call("rbind", lapply(comb, FUN = function ( y ) { ### 'dp' = statistics for derived parameters,siehe BIFIE-Hilfeseite von BIFIE.by
                           eval(parse(text=paste("dp <- list ( \"groupDiff\" =~ 0 + I(",y[1],"-",y[2],"))")))
                           resMd<- BIFIE.derivedParameters( resM, derived.parameters=dp )                 
   ### Achtung: falls 'group.differences.by' definiert, wird bereits auf der innersten Ebene, also quasi jetzt, begonnen, das wieder in die Ergebnisstruktur zurueckzupressen (ist, glaube ich, einfacher ... haha)
                           if ( length ( res ) == 1) {                          ### oh mann, das wird jetzt alles richtig bitter!
                                rg  <- "all.group=1"                            ### 'rg' = residual groups, hier wird unterschieden, ob es nur eine Gruppierungsvariable gibt 
                           }  else  {                                           ### (naemlich die, fuer die Gruppendifferenzen ausgegeben werden sollen), oder ob es mehrere
                                rg  <- setdiff(removeNonNumeric(res),0)         ### Gruppierungsvariablen gibt 
                                nam <- lapply(rg, FUN = function ( z ) { 
                                       nam <- liste[1,grep(paste0("groupvar", z), colnames(liste))]
                                       val <- x[1,paste0("groupval", z)]
                                       ret <- data.frame ( nam = nam, wert = val, toRec = dat.g[["labels"]][intersect(which(dat.g[["labels"]][,"varName"] == nam), which(dat.g[["labels"]][,"value"] == val)),"valLabel"], stringsAsFactors = FALSE)
                                       ret <- paste0(ret[["nam"]],"=",ret[["toRec"]])
                                       return(ret)})
                                rg  <- paste(unlist(nam), collapse=", ")
                           }                                                    ### Kontraste ausgeben; 'es' = effektstaerke
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
                           }                                                    ### untere Zeile: "est" und "es" werden mit minus 1 multipliziert, damit sie konsistent zu den survey-Ergebnissen sind
                           ret  <- data.frame ( group = paste(rg, vs, sep="___"), depVar = allNam[["dependent"]], modus = modus, comparison = "groupDiff", parameter = "mean", coefficient = c("est", "se", "p", "es"), value = c( (-1) * resMd[["stat"]][["coef"]],resMd[["stat"]][["se"]],resMd[["stat"]][["p"]], (-1)*es),add, stringsAsFactors = FALSE)
                           return(ret)}))
                   return(diffs)}))
      }  else  {
           grp  <- NULL
      }
   ### jetzt den Output in die Ergebnisstruktur zurueckpressen, Schritt 1: Namen der Gruppenvariablen umbenennen
      if(length(allNam[["group"]])==1) {                                        ### Boah, immer diese uneinheitlichen Spaltenbenennungen
           cols <- grep("groupva", colnames(resM[["stat"]]), value=TRUE)
           stopifnot(length(cols) ==2)
           altn <- data.frame ( alt = cols, neu = paste0(cols, "1"), stringsAsFactors = FALSE)
           recs <- paste("'",altn[,"alt"] , "' = '" , altn[,"neu"],"'",sep="", collapse="; ")
           colnames(resM[["stat"]]) <- recode(colnames(resM[["stat"]]), recs)
      }
      if(length(allNam[["group"]])>0) {
           for ( i in 1:length(allNam[["group"]])) {
                 labsG<- dat.g[["labels"]][which(dat.g[["labels"]][,"varName"] == allNam[["group"]][i]),]
                 recs <- paste("'",labsG[,"value"] , "' = '" , labsG[,"valLabel"],"'",sep="", collapse="; ")
                 resM[["stat"]][,paste0("groupval", i)] <- recode (resM[["stat"]][,paste0("groupval", i)], recs)
           }
      }
      if(!is.na(labsD[1,"value"])) {                                            ### nur fuer table: wenn die AV nicht numerisch war, musste sie zuvor fuer BIFIE numerisch gemacht werden
           recs <- paste("'",labsD[,"value"] , "' = '" , labsD[,"valLabel"],"'",sep="", collapse="; ")
           resM[["stat"]][,"varval"] <- recode(resM[["stat"]][,"varval"], recs) ### jetzt also wieder zuruecktransformieren in die urspruenglichen Werte
      }
      idv  <- c(grep("varval", colnames(resM[["stat"]]), value=TRUE), grep("groupval", colnames(resM[["stat"]]), value=TRUE))
      resML<- facToChar(melt(data=resM[["stat"]], id.vars=idv, measure.vars=mv,  na.rm=TRUE))
      resML[setdiff(1:nrow(resML), grep("_", resML[,"variable"])),"variable"] <- paste(resML[setdiff(1:nrow(resML), grep("_", resML[,"variable"])),"variable"], "est", sep="_")
      resML<- data.frame ( resML, colsplit(string = resML[,"variable"], pattern="_", names = c("parameter", "coefficient")),stringsAsFactors = FALSE)
      if ( toCall == "table") {
           resML <- resML[which(resML[,"parameter"] == "perc"),]
           resML[,"parameter"] <- resML[,"varval"]
           resML[,"variable"]  <- resML[,"varval"] <- NULL                      ### aus irgendeinem Grund darf in der Ergebnisstruktur keine Spalte stehen, die 'variable' heisst ... der Himmel weiss warum 
      }  else {
           resML[,"parameter"] <- recode(resML[,"parameter"], "'M'='mean'; 'SD'='sd'; 'Nweight'='NcasesValid'; 'Ncases'='sampleSize'")
           resML[,"variable"]  <- NULL
      }
      resML[,"coefficient"] <- recode(resML[,"coefficient"], "'SE'='se'")
      recs <- paste("'",grep("groupval", colnames(resML), value=TRUE) , "' = '" , allNam[["group"]],"'",sep="", collapse="; ")
      colnames(resML) <- recode(colnames(resML), recs)
      resML[,"modus"] <- modus
      resML[,"depVar"]<- allNam[["dependent"]]
      resML[,"comparison"] <- NA
      if ( length(allNam[["group"]]) > 1) {
           resML[,"group"] <- apply(resML[,allNam[["group"]]], MARGIN = 1, FUN = function ( z ) {paste(z, collapse="_")})
      }
      if ( length(allNam[["group"]]) == 1) { resML[,"group"] <- resML[,allNam[["group"]]]}
      resML<- rbind(resML, grp)
      return(resML)}

#bifie.table <- function( bifiObj , allNam, na.rm, group.delimiter, type, separate.missing.indicator, expected.values, modus){browser()}

#bifie.quantile <- function(bifiObj , allNam, na.rm, group.delimiter, type, probs, modus){browser()}

#bifie.glm <- function( bifiObj , allNam, formula, forceSingularityTreatment, glmTransformation, na.rm, group.delimiter, type, modus){browser()}

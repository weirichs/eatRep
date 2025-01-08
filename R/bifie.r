doBifieAnalyses <- function (dat.i, a){
      dat.i<- eatTools::facToChar(dat.i[,intersect(unlist(a[a%$$%allNam]), colnames(dat.i))], from = "character", to = "factor")
      dat.g<- eatGADS::import_DF(dat.i, checkVarNames = FALSE)                  ### obere zeile: aus performanzgruenden nur diejenigen variablen aus dem Datensatz auswaehlen, die spaeter auch gebraucht werden
      dat2 <- eatGADS::extractData(dat.g, convertLabels = "numeric")            ### Variablen muessen numerisch sein, geschieht hier ueber kleinen eatGADS Umweg
      dat2 <- sortBifie(dat = dat2, toCall=a$toCall, allNam=a[a%$$%allNam])
      labsD<- dat.g[["labels"]][which(dat.g[["labels"]][,"varName"] == a%$$%dependent),]
   ### in Listenformat bringen, to do: genestete Imputationen?
      datL <- by ( data = dat2, INDICES = dat2[,a%$$%imp], FUN = function ( imp.i ) { return(imp.i)})
      jkt  <- car::recode(a$type, "'JK2'='JK_TIMSS'; 'JK1'='JK_GROUP'")
   ### defaults fuer fayfac setzen. Fuer "JK_TIMSS" soll der Wert auf 2 gesetzt werden per default, das passiert aber innerhalb von BIFIE.data.jack, insofern wird der hier auf NULL belassen
	    if(is.null(a$rho)) {if(jkt == "JK_GROUP") {a$rho <- (length(unique(dat.i[, a%$$%PSU]))-1)/length(unique(dat.i[, a%$$%PSU]))   }}
      txt  <- capture.output(bo   <- BIFIEsurvey::BIFIE.data.jack( data= datL,  wgt = a%$$%wgt, jktype=jkt , jkzone = a%$$%PSU, jkrep = a%$$%repInd, jkfac=a%$$%jkfac, fayfac=a%$$%rho, cdata=FALSE))
      if ( isTRUE(a%$$%verbose)) { cat("\n"); printJackInfo(bo=bo, a=a, jkt=jkt, cdata=FALSE)}
      attributes(a$group) <- NULL                                               ### Attribute der Gruppierungsvariablen entfernen, sonst gibt BIFIEsurvey einen Fehler aus
   ### jetzt analysespezifische Subfunktion starten
      if ( a%$$%toCall == "mean") {
           resML <- bifieMean(bifie.obj = bo, allNam=a[a%$$%allNam], dat.g=dat.g, labsD=labsD, toCall=a%$$%toCall, dat.i=dat.i, modus=a%$$%modus, verbose=a%$$%verbose)
      }
      if ( a%$$%toCall == "table") {
           resML <- bifieTable(bifie.obj = bo, allNam=a[a%$$%allNam], dat.g=dat.g, labsD=labsD, toCall=a%$$%toCall, dat.i=dat.i, modus=a%$$%modus, verbose=a%$$%verbose)
      }
      if ( a%$$%toCall == "lmer") {
           foo   <- checkWithinBetweenWeights(dat=dat.i, allNam=a[a%$$%allNam])
           resML <- bifieLmer(bifie.obj=bo, allNam=a[a%$$%allNam], dat.g=dat.g, labsD=labsD, modus=a%$$%modus, formula.fixed=a%$$%formula.fixed, formula.random=a%$$%formula.random, verbose=a%$$%verbose)
      }
      return(resML)}

sortBifie <- function(dat, toCall, allNam) {
      if ( toCall != "lmer") {
           dat <- dat[order(dat[,allNam[["ID"]]]), ]                            ### fuer alles andere als multilevel analysen: sortieren nach ID
      }  else  {
           dat <- dat[order(dat[,allNam[["clusters"]]]), ]                      ### fuer multilevel analysen: sortieren nach clustern
      }
      return(dat)}

### bifie.survey output in die Ergebnisstruktur zuruecktransformieren
aufbOut <- function (resM, allNam, dat.g, labsD) {
   ### Schritt 1: Namen der Gruppenvariablen umbenennen
      if(length(allNam[["group"]])==1) {
           cols <- grep("groupva", colnames(resM[["stat"]]), value=TRUE)
           stopifnot(length(cols) ==2)
           altn <- data.frame ( alt = cols, neu = paste0(cols, "1"), stringsAsFactors = FALSE)
           colnames(resM[["stat"]]) <- eatTools::recodeLookup(colnames(resM[["stat"]]), altn)
      }
      if(length(allNam[["group"]])>0) {
           for ( i in 1:length(allNam[["group"]])) {
                 labsG<- dat.g[["labels"]][which(dat.g[["labels"]][,"varName"] == allNam[["group"]][i]),]
                 resM[["stat"]][,paste0("groupval", i)] <- eatTools::recodeLookup(resM[["stat"]][,paste0("groupval", i)], labsG[,c("value", "valLabel")])
           }
      }
      if(!is.na(labsD[1,"value"])) {                                            ### nur fuer table: wenn die AV nicht numerisch war, musste sie zuvor fuer BIFIE numerisch gemacht werden
           recs <- paste("'",labsD[,"value"] , "' = '" , labsD[,"valLabel"],"'",sep="", collapse="; ")
           resM[["stat"]][,"varval"] <- car::recode(resM[["stat"]][,"varval"], recs)
      }                                                                         ### jetzt also wieder zuruecktransformieren in die urspruenglichen Werte
      return(resM)}

### Aufbereitung fuer NICHT-multilevel-analysen
aufbNoMultilevel <- function (resM, mv, toCall, allNam, dat.i, modus, grp) {
          idv  <- c(grep("varval", colnames(resM[["stat"]]), value=TRUE), grep("groupval", colnames(resM[["stat"]]), value=TRUE))
          resML<- eatTools::facToChar(reshape2::melt(data=resM[["stat"]], id.vars=idv, measure.vars=mv,  na.rm=TRUE))
          resML[setdiff(1:nrow(resML), grep("_", resML[,"variable"])),"variable"] <- paste(resML[setdiff(1:nrow(resML), grep("_", resML[,"variable"])),"variable"], "est", sep="_")
          resML<- data.frame ( resML, reshape2::colsplit(string = resML[,"variable"], pattern="_", names = c("parameter", "coefficient")),stringsAsFactors = FALSE)
          if ( toCall == "table") {
               resML <- resML[which(resML[,"parameter"] == "perc"),]
               resML[,"parameter"] <- resML[,"varval"]
               resML[,"variable"]  <- resML[,"varval"] <- NULL                  ### in der Ergebnisstruktur darf keine Spalte stehen, die 'variable' heisst
          }  else {
               resML[,"parameter"] <- car::recode(resML[,"parameter"], "'M'='mean'; 'SD'='sd'; 'Nweight'='NcasesValid'")
               resML[,"variable"]  <- NULL
          }
          resML[,"coefficient"] <- car::recode(resML[,"coefficient"], "'SE'='se'")
          recs <- paste("'",grep("groupval", colnames(resML), value=TRUE) , "' = '" , allNam[["group"]],"'",sep="", collapse="; ")
          colnames(resML) <- car::recode(colnames(resML), recs)
          if ( toCall == "table") {                                             ### Stichprobengroessen fuer table dazu
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

bifieMean <- function ( bifie.obj, allNam, dat.g, labsD, toCall, dat.i, modus, verbose) {
      if(verbose) {cat(paste0("'BIFIE.univar' for 'call = mean'. dependent = '",allNam[["dependent"]],"'. group(s) = '",paste(allNam[["group"]], collapse="', '"), "'. \n"))}
      txt  <- capture.output(resM <- BIFIEsurvey::BIFIE.univar( BIFIEobj=bifie.obj , vars = allNam[["dependent"]], group=allNam[["group"]] ))
      mv   <- c("Nweight", "Ncases", "M", "M_SE", "M_p", "SD", "SD_SE", "SD_p")
      grp  <- computeGroupDifferences(resM=resM, allNam=allNam, dat.g=dat.g, modus=modus)
      resM <- aufbOut(resM=resM, allNam=allNam, dat.g=dat.g, labsD=labsD)
      resML<- aufbNoMultilevel(resM=resM, mv=mv, toCall=toCall, allNam=allNam, dat.i=dat.i, modus=modus, grp=grp)
      return(resML)  }

bifieTable <- function ( bifie.obj, allNam, dat.g, labsD, toCall, dat.i, modus, verbose) {
      if(verbose) {cat(paste0("'BIFIE.freq' for 'call = table'. dependent = '",allNam[["dependent"]],"'. group(s) = '",paste(allNam[["group"]], collapse="', '"), "'. \n"))}
      txt  <- capture.output(resM <- BIFIEsurvey::BIFIE.freq( BIFIEobj=bifie.obj , vars = allNam[["dependent"]], group=allNam[["group"]] ))
      resM[["stat"]][,"perc_p"] <- 2*(1-pnorm(resM[["stat"]][,"perc"] / resM[["stat"]][,"perc_SE"]))
      mv   <- c("Nweight", "Ncases", "perc", "perc_SE", "perc_p")               ### measure.vars fuers reshapen
      stopifnot(is.null(allNam[["group.differences.by"]]))                      ### fuer table darfs keine group.differences.by geben
      resM <- aufbOut(resM=resM, allNam=allNam, dat.g=dat.g, labsD=labsD)
      resML<- aufbNoMultilevel(resM=resM, mv=mv, toCall=toCall, allNam=allNam, dat.i=dat.i, modus=modus, grp=NULL)
      return(resML)  }
      
bifieLmer <- function ( bifie.obj, allNam, dat.g, labsD, modus, formula.fixed, formula.random, verbose) {
      resM <- BIFIEsurvey::BIFIE.twolevelreg( BIFIEobj=bifie.obj, dep=allNam[["dependent"]], formula.fixed=formula.fixed, formula.random=formula.random, idcluster=allNam[["clusters"]], wgtlevel2=allNam[["L2wgt"]], wgtlevel1=allNam[["L1wgt"]], group = allNam[["group"]])
      resM <- aufbOut(resM=resM, allNam=allNam, dat.g=dat.g, labsD=labsD)       ### untere Zeile: Aufbereitung fuer multilevel analysen
      resML<- aufbMultilevel(resM=resM, allNam=allNam, modus=modus)
      return(resML)}

### nur fuer die mean-funktion
computeGroupDifferences <- function(resM, allNam, dat.g, modus){
      if(is.null(allNam[["group.differences.by"]])) {return(NULL)} else {       ### die 'parnames' haben nur fuer die BIFIE.univar-Funktion komische krude Namen, deswegen ist das hier alles komplizierter als es die Hilfeseite von BIFIE.by andeutet ... und man muss UNBEDINGT mit sort=FALSE mergen, sonst stimmt die Reihenfolge in 'resM' nicht mit der in 'redMd' ueberein!
           liste<- data.frame ( merge(resM[["stat_M"]],resM[["stat_SD"]][,c(grep("^groupva", colnames(resM[["stat_SD"]]), value=TRUE), "SD", "SD_SE")],  by = grep("^groupva", colnames(resM[["stat_SD"]]), value=TRUE), all=TRUE, sort=FALSE), dp = resM[["parnames"]], groupvar0 = "wholePop", groupval0 = 0, stringsAsFactors = FALSE)
           if ( length(allNam[["group.differences.by"]]) == length(allNam[["group"]])) {
               colnames(liste) <- car::recode(colnames(liste), "'groupvar'='groupvar1'; 'groupval'='groupval1'")
           }
           col  <- colnames(liste)[which(sapply(eatTools::facToChar(liste[1,]), FUN = function ( x ) { x == allNam[["group.differences.by"]] }))]
           col  <- c(eatTools::removeNumeric(col), eatTools::removeNonNumeric(col))
           res  <- setdiff(grep("^groupval", colnames(liste), value=TRUE), paste0("groupval", col[2]))
           grp  <- do.call("rbind", by(data=liste, INDICES = liste[,res], FUN = function ( x ) {
                   comb <- data.frame ( combinat::combn(x=x[,"dp"], m=2), stringsAsFactors = FALSE)
                   diffs<- do.call("rbind", lapply(comb, FUN = function ( y ) { ### 'dp' = statistics for derived parameters,siehe BIFIE-Hilfeseite von BIFIE.by
                           dp   <- eval(parse(text=paste("list ( \"groupDiff\" =~ 0 + I(",y[1],"-",y[2],"))")))
                           resMd<- BIFIEsurvey::BIFIE.derivedParameters( resM, derived.parameters=dp )
   ### Achtung: falls 'group.differences.by' definiert, wird bereits auf der innersten Ebene, also quasi jetzt, begonnen, das wieder in die Ergebnisstruktur zurueck zu ueberfuehren
                           if ( length ( res ) == 1) {
                                rg  <- "all.group=1"                            ### 'rg' = residual groups, hier wird unterschieden, ob es nur eine Gruppierungsvariable gibt
                           }  else  {                                           ### (naemlich die, fuer die Gruppendifferenzen ausgegeben werden sollen), oder ob es mehrere
                                rg  <- setdiff(eatTools::removeNonNumeric(res),0)## Gruppierungsvariablen gibt
                                nam <- lapply(rg, FUN = function ( z ) {
                                       nam <- liste[1,grep(paste0("groupvar", z), colnames(liste))]
                                       val <- x[1,paste0("groupval", z)]
                                       ret <- data.frame ( nam = nam, wert = val, toRec = dat.g[["labels"]][intersect(which(dat.g[["labels"]][,"varName"] == nam), which(dat.g[["labels"]][,"value"] == val)),"valLabel"], stringsAsFactors = FALSE)
                                       ret <- paste0(ret[["nam"]],"=",ret[["toRec"]])
                                       return(ret)})
                                rg  <- paste(unlist(nam), collapse=", ")
                           }                                                    ### Kontraste ausgeben; 'es' = effektstaerke
                           vs   <- x[which(x[,"dp"] %in% y),paste(c("groupval",col[2]),collapse="")]
                           vs   <- paste(dat.g[["labels"]][intersect(which(dat.g[["labels"]][,"varName"] %in% allNam[["group.differences.by"]]), which(dat.g[["labels"]][,"value"] %in% vs)),"valLabel"], collapse=" - ")
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
                           ret  <- data.frame ( group = paste(rg, vs, sep="___"), depVar = allNam[["dependent"]], modus = paste(modus,"BIFIEsurvey", sep="__"),  comparison = "groupDiff", parameter = "mean", coefficient = c("est", "se", "p", "es"), value = c( (-1) * resMd[["stat"]][["coef"]],resMd[["stat"]][["se"]],resMd[["stat"]][["p"]], (-1)*es),add, stringsAsFactors = FALSE)
                           return(ret)}))
                   return(diffs)}))
      }
      return(grp)}

### Aufbereitung fuer multilevel analysen
aufbMultilevel <- function(resM, allNam, modus) {
      if(!is.null(allNam[["group"]])) {
          colnames(resM[["stat"]]) <- car::recode(colnames(resM[["stat"]]), paste0("'groupval1'='",allNam[["group"]],"'"))
          resM[["stat"]][,"group"] <- resM[["stat"]][,allNam[["group"]]]
      }  else  {
          resM[["stat"]][,"group"] <- "wholeGroup"
      }                                                                         ### jetzt reshapen und dann weiter getrennt fuer mit/ohne Gruppenvariable aufbereiten
      resL <- eatTools::facToChar(reshape2::melt(data=resM[["stat"]], id.vars=c("parameter", "group", allNam[["group"]]), measure.vars=c("est", "SE", "p"),  na.rm=TRUE))
   ### to do: hier lieber alle Parameter, also auch random slope varianz etc. so mitnehmen wie sie sind und nur fuer die trends ausschliesslich die regressionskoeffizienten auswaehlen!
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

checkWithinBetweenWeights <- function(dat, allNam){
   ### Anzahl an Gesamtgewichten sollte nicht geringer sein als L1- oder L2-Gewichte
      if(!is.null(allNam[["L2wgt"]]) && !is.null(allNam[["wgt"]])) {
         if (!is.null(allNam[["imp"]])) {
             d <- dat[which(dat[,allNam[["imp"]]] == sort(dat[,allNam[["imp"]]])[1]),]
         }  else  {
             d <- dat
         }
         nUnit <- sapply( d[,c(allNam[["L1wgt"]], allNam[["L2wgt"]], allNam[["wgt"]])], FUN = function (x) {length(unique(x))})
         if (!is.null(allNam[["L1wgt"]])) {if ( nUnit[[allNam[["wgt"]]]] < nUnit[[allNam[["L1wgt"]]]]) {warning(paste0("Number of distinct units in sampling weights '",allNam[["wgt"]],"' (",nUnit[[allNam[["wgt"]]]],"), is lower than the number of distinct units in level-1 weights '",allNam[["L1wgt"]],"' (",nUnit[[allNam[["L1wgt"]]]],")"))}}
         if (!is.null(allNam[["L2wgt"]])) {if ( nUnit[[allNam[["wgt"]]]] < nUnit[[allNam[["L2wgt"]]]]) {warning(paste0("Number of distinct units in sampling weights '",allNam[["wgt"]],"' (",nUnit[[allNam[["wgt"]]]],"), is lower than the number of distinct units in level-2 weights '",allNam[["L2wgt"]],"' (",nUnit[[allNam[["L2wgt"]]]],")"))}}
      }}

printJackInfo <- function(bo, a, jkt, cdata){
      pre  <- match.call(BIFIEsurvey::BIFIE.data.jack, call("BIFIEsurvey::BIFIE.data.jack", "datL", wgt = a%$$%wgt, jktype=jkt , jkzone = a%$$%PSU, jkrep = a%$$%repInd, jkfac=a%$$%jkfac, fayfac=a%$$%rho, cdata=FALSE))
      txt1 <- capture.output(print(bo))
      start<- min(grep("data with", txt1))
      print(pre)
      cat(txt1[5:length(txt1)], sep = " || "); cat("\n")}

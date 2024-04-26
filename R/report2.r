### setzt auf dem Output von report() auf
report2 <- function(repFunOut, trendDiffs = FALSE, add=list(), exclude = c("NcasesValid", "var", "sampleSize"), printGlm = FALSE,
                     round = TRUE, digits = 3, printDeviance = FALSE) {
           out   <- eval(parse(text=paste0("report(",paste(formalArgs(report), formalArgs(report), collapse=", ", sep="="), ")")))
           measur<- setdiff(colnames(out), c(intersect(unlist(attr(out, "allNam")), colnames(out)), names(add), c("group", "depVar", "modus", "comparison", "parameter")))
           outL  <- reshape2::melt(out, measure.vars = measur, na.rm=TRUE)
           if(length(grep("_", as.character(outL[,"variable"])))>0) {
              outL[,"coef"] <- car::recode(eatTools::halveString(as.character(outL[,"variable"]), "_", first=TRUE)[,1], "'sig'='p'")
           }  else  {
              outL[,"coef"] <- as.character(outL[,"variable"])
           }
           if(!is.null(attr(out, "allNam")[["trend"]])) {
              outL[,attr(out, "allNam")[["trend"]]] <- eatTools::crop(eatTools::removePattern(eatTools::halveString(as.character(outL[,"variable"]), "_", first=TRUE)[,2], "trend"), "_")
           }
           outL[,"variable"] <- NULL
           if(any(c("trendDiff_cross", "trendDiff_group") %in% outL[,"comparison"])) {
              outL  <- outL[-eatTools::whereAre(c("trendDiff_cross", "trendDiff_group"), outL[,"comparison"], verbose=FALSE),]
           }
           if(!is.null(attr(out, "allNam")[["trend"]])) {
              isTren<- grep(".vs.", outL[,attr(out, "allNam")[["trend"]]])
              outL[isTren,"comparison"] <- eatTools::crop(paste("trend", car::recode(outL[isTren,"comparison"],"NA=''"), sep="_"), "_")
           }                                                                    ### untere Zeile: das hier soll spaeter im tabellenblatt 'plain' landen
           plain <- eatTools::makeDataFrame(tidyr::pivot_wider(outL, names_from = "coef", values_from = "value"))
           plain[,"row"] <- 1:nrow(plain)
    ### Schritt 1: id fuer comparison = NA vergeben
           plainn<- subset(plain, is.na(comparison))
           plainn<- do.call("rbind", by(data = plainn, INDICES = plainn[,c("group", "parameter", attr(out, "allNam")[["trend"]], "modus", "depVar")], FUN = function (sg) {
                    stopifnot(length(sg[,"parameter"]) == length(unique(sg[,"parameter"])))
                    sg[,"id"] <- paste("group", min(sg[,"row"]), sep="_")
                    return(sg)}))
           plainc<- subset(plain, !is.na(comparison))
           comps <- gsub("_of_", "_", unique(plainc[,"comparison"]))
           comps <- data.frame ( comp = comps, compOri = unique(plainc[,"comparison"]), level = sapply(strsplit(comps, "_"), length), stringsAsFactors = FALSE)
           if(nrow(plainc)>0 && nrow(subset(comps, level == 1))>0) {
               comps1<- do.call("rbind", apply(subset(comps, level == 1), MARGIN = 1, FUN = function (co1) {
                        plc1 <- subset(plainc, comparison == co1[["comp"]])         ### untere Zeile: jetzt die entsprechenden Vergleichseintraege finden!
                        plc1i<- do.call("rbind", by(data = plc1, INDICES = plc1[,c("group", "parameter", attr(out, "allNam")[["trend"]], "modus", "depVar")], FUN = function (sg1) {
                                if ( unique(sg1[,"parameter"]) == "chiSquareTest") {return(sg1)}
                                fun <- paste0("idFocRefGrp",toupper(sg1[,"comparison"]))
                                foc <- eval(parse(text=paste0(fun, "(string = sg1, out=out, plainn=plainn, add=add)")))
                                sg1 <- data.frame ( sg1, id = paste("comp", sg1[,"row"], sep="_"), unit_1 = foc[1,"id"], unit_2 = foc[2,"id"], stringsAsFactors = FALSE)
                                return(sg1)}))
                        return(plc1i)}, simplify = FALSE))
           }  else  {
               comps1<- NULL
           }
    ### Schritt 2: id fuer Zweifach-Gruppenvergleiche vergeben
           if(nrow(plainc)>0 && nrow(subset(comps, level == 2))>0) {
               comps2<- do.call("rbind", apply(subset(comps, level == 2), MARGIN = 1, FUN = function (co1) {
                        plc2 <- subset(plainc, comparison == co1[["compOri"]])
                        plc2i<- do.call("rbind", by(data = plc2, INDICES = plc2[,c("group", "parameter", attr(out, "allNam")[["trend"]], "modus", "depVar")], FUN = function (sg2) {
                                fun <- paste0("idFocRefGrp",car::recode(toupper(sg2[,"comparison"]), "'TREND_CROSSDIFF' = 'TREND'; 'TREND_GROUPDIFF'='TREND'"))
                                foc <- eval(parse(text=paste0(fun, "(string = sg2, out=out, plainn=comps1, add=add)")))
                                sg2 <- data.frame ( sg2, id = paste("comp", sg2[,"row"], sep="_"), unit_1 = foc[1,"id"], unit_2 = foc[2,"id"], stringsAsFactors = FALSE)
                                return(sg2)}))
                        return(plc2i)}, simplify = FALSE))
           }  else  {
               comps2<- NULL
           }
    ### Schritt 3: id fuer drefach-Gruppenvergleiche vergeben
           if(nrow(plainc)>0 && nrow(subset(comps, level == 3))>0) {
               comps3<- do.call("rbind", apply(subset(comps, level == 3), MARGIN = 1, FUN = function (co1) {
                        plc3 <- subset(plainc, comparison == co1[["compOri"]])
                        plc3i<- do.call("rbind", by(data = plc3, INDICES = plc3[,c("group", "parameter", attr(out, "allNam")[["trend"]], "modus", "depVar")], FUN = function (sg3) {
                                fun <- paste0("idFocRefGrp",car::recode(toupper(sg3[,"comparison"]), "'TREND_CROSSDIFF' = 'TREND'; 'TREND_GROUPDIFF'='TREND'; 'TREND_CROSSDIFF_OF_GROUPDIFF'='TREND'"))
                                foc <- eval(parse(text=paste0(fun, "(string = sg3, out=out, plainn=comps2, add=add)")))
                                sg3 <- data.frame ( sg3, id = paste("comp", sg3[,"row"], sep="_"), unit_1 = foc[1,"id"], unit_2 = foc[2,"id"], stringsAsFactors = FALSE)
                                return(sg3)}))
                        return(plc3i)}, simplify = FALSE))
           }  else  {
               comps3<- NULL
           }
    ### Schritt 4: alles zusammenbaemmeln
           plain <- plyr::rbind.fill(plainn, comps1, comps2, comps3)
           miss  <- which(is.na(plain[,"id"]))
           if(length(miss)>0) {
               type <- unique(plain[miss,"parameter"])
               stopifnot(length(setdiff("chiSquareTest", type))==0)
           }
    ### Schritt 5: die sheets bauen entsprechend 'p:\Methoden\02_R_Pakete\eatRep\Outputstruktur\eatRep_demo.xlsx'
           compar<- subset(plain, !is.na(comparison))[,na.omit(match(c("id", "unit_1", "unit_2", "comparison"), colnames(plain)))]
           groups<- subset(plain, is.na(comparison))[,na.omit(match(c("id", attr(out, "allNam")[["group"]], attr(out, "allNam")[["trend"]],  names(add)), colnames(plain))), drop=FALSE]
           if(length(attr(out, "allNam")[["group"]])>0) {for ( i in attr(out, "allNam")[["group"]]) {if(i %in% colnames(groups)) {groups[,i] <- car::recode(groups[,i], "NA='all'")}}}
           estim <- plain[,c("id", colnames(plain)[na.omit(match(c("depVar", "parameter", "est", "se", "p", "es"), colnames(plain)))])]
    ### Schritt 6: Anpassungen entsprechend der Mail von Nicklas, 21.03.2024, 8.32 Uhr: "gepasteten strings" anpassen
           plain1<- rbind(subset(plain, is.na(comparison)), subset(plain, comparison=="trend"))
           if(length(attr(out, "allNam")[["group"]]) >0) {                      ### fuer alle Levels der comparison variable getrennt
              plain1[,"group"] <- apply(plain1[,attr(out, "allNam")[["group"]], drop=FALSE], MARGIN = 1, FUN = function (zeile) {
                                  z <- na.omit(zeile)
                                  return(car::recode(paste(names(z), z, sep="=", collapse=", "), "''='total'"))})
           }
           plain2<- subset(plain, comparison %in% c("trend_crossDiff", "crossDiff"))
           if(nrow(plain2)>0) {
              plain2 <- do.call("rbind", by(data = plain2, INDICES = plain2[,"group"], FUN = function (grp) {
                        halb <- strsplit(unique(grp[,"group"]), ".vs.")[[1]]
                        halb1<- lapply(halb, FUN = function (x) {
                                str1 <- strsplit(x, split="_")[[1]]
                                str1 <- addVarNamsToVarVals (string=str1, group = attr(out, "allNam")[["group"]], plain = plain1, tot = "total") })
                        grp[,"group"]<- paste(paste(unlist(halb1[[1]]), collapse= ", ") , paste(unlist(halb1[[2]]), collapse= ", "), sep = " - ")
                        for ( i in 1:length(halb1)) {for ( j in 1:length(halb1[[i]])) {if(!is.null(attr(halb1[[i]][[j]], "dfr"))) {if(is.na(grp[1,attr(halb1[[i]][[j]], "dfr")[["variable"]]])) {grp[,attr(halb1[[i]][[j]], "dfr")[["variable"]]] <- paste0(attr(halb1[[i]][[j]], "dfr")[["value"]], ".vs.total") }}}}
                        return(grp)}))
           }
           plain3<- subset(plain, comparison %in% c("trend_groupDiff", "groupDiff"))
           if(nrow(plain3)>0) {
              second <- eatTools::halveString(plain3[,"group"], "____")         ### untere Zeile: Achtung! hier muss wirklich 'plain1' stehen, weil 'plain3' keine Gruppenvalues fuer crossdiff enthaelt
              txt    <- capture.output(third <- eatTools::makeDataFrame(eatTools::halveString(second[,2], ".vs.")))
              third  <- data.frame ( lapply(third, FUN = function (col) {unlist(addVarNamsToVarVals (string=col, group = attr(out, "allNam")[["group"]], plain = plain1, tot = "total"))}), stringsAsFactors = FALSE)
              plain3[,"group"] <- paste0(gsub("all.group=1", "total", second[,1]), ": ", third[,1], " - ", third[,2])
           }
           plain4<- subset(plain, comparison %in% c("crossDiff_of_groupDiff", "trend_crossDiff_of_groupDiff"))
           if(nrow(plain4)>0) {
              plain4 <- do.call("rbind", by(data = plain4, INDICES = plain4[,"group"], FUN = function (grp) {
                        spl1   <- strsplit(unique(grp[,"group"]), ".VS.")[[1]]
                        spl2   <- strsplit(spl1, "____")
                        allVals<- unlist(strsplit(unlist(spl2), "=|.vs."))
                        allStrg<- addVarNamsToVarVals (string=allVals, group = attr(out, "allNam")[["group"]], plain = plain1, tot = "total")
                        for ( i in 1:length(allStrg)) {if(!is.null(attr(allStrg[[i]], "dfr"))) {if(is.na(grp[1,attr(allStrg[[i]], "dfr")[["variable"]]])) {grp[,attr(allStrg[[i]], "dfr")[["variable"]]] <- paste0(attr(allStrg[[i]], "dfr")[["value"]], ".vs.total") }}}
                        spl2   <- lapply(spl2, FUN = function (x) {
                                  tx<- capture.output(y <- eatTools::makeDataFrame(eatTools::halveString(x[2], ".vs.")))
                                  y <- lapply(y, FUN = function (col) {addVarNamsToVarVals (string=col, group = attr(out, "allNam")[["group"]], plain = plain1, tot = "total")})
                                  y <- gsub("all.group=1", "total", c(x[1], unlist(y[[1]]), unlist(y[[2]])))
                                  stopifnot(length(y) == 3)
                                  y <- paste0("(",y[1], ": ",y[2], " - ", y[3], ")")
                                  return(y)})
                        stopifnot(length(spl2) ==2)
                        grp[,"group"] <- paste0(spl2[1], " - ", spl2[2])
                        return(grp)}))
           }
           stopifnot(nrow(plain) == nrow(rbind(plain1, plain2, plain3, plain4)))
           plain <- rbind(plain1, plain2, plain3, plain4)
           rownames(plain) <- NULL
           return(list(plain = plain[,-match("row", colnames(plain))], comparisons = compar, group=groups, estimate = estim))}

### Funktion waehlt die Fokus- und Referenzgruppe fuer cross differences aus
idFocRefGrpCROSSDIFF <- function(string, out, plainn, add) {
           prepost<- strsplit(unique(string[,"group"]), ".vs.")
           stopifnot(length(prepost[[1]]) ==2)
           prepost<- lapply(prepost, FUN = function (x) {strsplit(x, "_")})
           focCols<- data.frame ( lapply(plainn, FUN = function (col) { any(prepost[[1]][[1]] %in% col) })[attr(out, "allNam")[["group"]]], drop=FALSE, stringsAsFactors = FALSE)
           focCols<- focCols[,which(sapply(focCols, FUN = eval)), drop=FALSE]
           for ( i in 1:ncol(focCols)) {focCols[,i] <- prepost[[1]][[1]][i]}
           focCols<- data.frame ( focCols, string[,c(attr(out, "allNam")[["trend"]], "parameter", "depVar")],stringsAsFactors = FALSE)
           focus  <- merge(focCols,plainn, by = colnames(focCols), all=FALSE)
           if(nrow(focus) != 1) {focus  <- chooseLine(obj = focus, prepost=prepost, whichOne = 1, char = "_", groupName = "group")}
           if( prepost[[1]][[2]][1] == "wholeGroup") {
               refCol <- data.frame ( group = "wholeGroup", string[,c(attr(out, "allNam")[["trend"]], "parameter", "depVar")], stringsAsFactors = FALSE)
               refere <- merge(refCol,plainn, by = colnames(refCol), all=FALSE)
           }  else  {
               refCol <- data.frame ( lapply(plainn, FUN = function (col) { any(prepost[[1]][[2]] %in% col) })[attr(out, "allNam")[["group"]]], drop=FALSE, stringsAsFactors = FALSE)
               refCol <- refCol[,which(sapply(refCol, FUN = eval)), drop=FALSE]
               for ( i in 1:ncol(refCol)) {refCol[,i] <- prepost[[1]][[2]][i]}
               refCol <- data.frame ( refCol, string[,c(attr(out, "allNam")[["trend"]], "parameter", "depVar")], stringsAsFactors = FALSE)
               refere <- merge(refCol,plainn, by = colnames(refCol), all=FALSE)
               if ( length(prepost[[1]][[2]]) == 1) {refere  <- refere[setdiff(1:nrow(refere), grep("_", refere[,"group"])),]}
           }
           if(nrow(refere) != 1) {refere <- chooseLine(obj = refere, prepost=prepost, whichOne = 2, char = "_", groupName = "group")}
           return(rbind(focus, refere))}
           
chooseLine <- function(obj, prepost, whichOne, char, groupName) {
           anzahl <- stringr::str_count(obj[,groupName], char)
           needed <- length(prepost[[1]][[whichOne]]) - 1
           stopifnot ( length(which(anzahl == needed)) == 1)
           obj    <- obj[which(anzahl == needed),]
           stopifnot(nrow(obj) ==1)
           return(obj)}

### Funktion waehlt die Fokus- und Referenzgruppe fuer trend differences aus
idFocRefGrpTREND <- function(string, out, plainn, add) {
           tomerge<- string[,c(attr(out, "allNam")[["trend"]], "group", "parameter", "modus", "depVar", names(add))]
           tomerge<- rbind(tomerge, tomerge)
           tomerge[,1] <- strsplit(tomerge[1,1], ".vs.")[[1]]
           tomerge<- merge(tomerge,plainn, by = colnames(tomerge), all=FALSE)
           stopifnot(nrow(tomerge)==2)
           return(tomerge)}

### Funktion waehlt die Fokus- und Referenzgruppe fuer group differences aus
idFocRefGrpGROUPDIFF <- function(string, out, plainn, add) {
           kontras<- unlist(lapply(string[,attr(out, "allNam")[["group"]], drop=FALSE], FUN = function (col) {grep(".vs.", col)}))
           stopifnot(length(kontras) == 1)
           kontr  <- strsplit(string[,names(kontras)], split = ".vs.")
           tomerge<- string[,c(attr(out, "allNam")[["trend"]], "parameter", "modus", "depVar", attr(out, "allNam")[["group"]], names(add))]
           tomerge<- rbind(tomerge, tomerge)
           tomerge[,names(kontras)] <- kontr[[1]]
           tomerge<- merge(tomerge,plainn, by = colnames(tomerge), all=FALSE)
           stopifnot(nrow(tomerge)==2)
           return(tomerge)}

idFocRefGrpCROSSDIFF_OF_GROUPDIFF <- function (string, out, plainn, add) {
    ### erstmal den Gruppenvergleich finden (Referenz!)
           spl1   <- strsplit(string[,"group"], split = ".VS.")
           tomerge<- string[,c(attr(out, "allNam")[["trend"]], "parameter", "modus", "depVar", names(add))]
           tomerge[,"group"] <- spl1[[1]][2]
           refere <- merge(tomerge, subset(plainn, comparison == "groupDiff"), by = colnames(tomerge), all=FALSE)
           stopifnot(nrow(refere)==1)
    ### jetzt den crossdiff finden
           focus  <- string[,c(attr(out, "allNam")[["trend"]], "parameter", "modus", "depVar", names(add))]
           focus[,"group"] <- spl1[[1]][1]
           focus  <- merge(focus, subset(plainn, comparison == "groupDiff"), by = colnames(focus), all=FALSE)
           stopifnot(nrow(refere)==1)
           return(rbind(focus, refere))}

addVarNamsToVarVals <- function (string, group , plain , tot = "total"){
           str1 <- lapply(string, FUN = function (st) {
                   col <- sapply(plain[,group, drop=FALSE], FUN = function (x) {st %in% x})
                   col <- names(which(col==TRUE))
                   if ( length(col) ==0) {return(tot)}
                   dfr <- list(variable = col, value = st)
                   ret <- paste0(col, "=", st)
                   attr(ret, "dfr") <- dfr
                   return(ret)})
           return(str1)}



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
           plainn<- do.call("rbind", by(data = plainn, INDICES = plainn[,c("group", attr(out, "allNam")[["trend"]], "depVar")], FUN = function (sg) {
                    stopifnot(length(sg[,"parameter"]) == length(unique(sg[,"parameter"])))
                    sg[,"id"] <- paste("group", min(sg[,"row"]), sep="_")
                    return(sg)}))
           plainc<- subset(plain, !is.na(comparison))
           comps <- gsub("_of_", "_", unique(plainc[,"comparison"]))
           comps <- data.frame ( comp = comps, compOri = unique(plainc[,"comparison"]), level = sapply(strsplit(comps, "_"), length), stringsAsFactors = FALSE)
           comps1<- setIDs(comps, lev=1, plainc=plainc,oldNew =NULL,plainn=plainn, co1 = "comp", indices = c("group", attr(out, "allNam")[["trend"]], "depVar"), out=out, add=add)
    ### Schritt 2: id fuer Zweifach-Gruppenvergleiche vergeben
           comps2<- setIDs(comps, lev=2, plainc=plainc,oldNew =data.frame ( old = c("TREND_CROSSDIFF", "TREND_GROUPDIFF"), new = rep("TREND", 2), stringsAsFactors = FALSE), plainn=comps1, co1 = "compOri", indices = c("group", attr(out, "allNam")[["trend"]], "depVar"), out=out, add=add)
    ### Schritt 3: id fuer dreifach-Gruppenvergleiche vergeben
           comps3<- setIDs(comps, lev=3, plainc=plainc,oldNew =data.frame ( old = c("TREND_CROSSDIFF", "TREND_GROUPDIFF", "TREND_CROSSDIFF_OF_GROUPDIFF"), new = rep("TREND", 3), stringsAsFactors = FALSE), plainn=comps2, co1 = "compOri", indices = c("group", "parameter", attr(out, "allNam")[["trend"]], "modus", "depVar"), out=out, add=add)
    ### Schritt 4: alles zusammenbaemmeln
           plain <- plyr::rbind.fill(plainn, comps1, comps2, comps3)
           miss  <- which(is.na(plain[,"id"]))
           if(length(miss)>0) {
               type <- unique(plain[miss,"parameter"])
               stopifnot(length(setdiff("chiSquareTest", type))==0)
           }
    ### Schritt 5: die sheets bauen entsprechend 'p:\Methoden\02_R_Pakete\eatRep\Outputstruktur\eatRep_demo.xlsx'
           compar<- unique(subset(plain, !is.na(comparison))[,na.omit(match(c("id", "unit_1", "unit_2", "comparison"), colnames(plain)))])
           groups<- unique(subset(plain, is.na(comparison))[,na.omit(match(c("id", attr(out, "allNam")[["group"]], attr(out, "allNam")[["trend"]],  names(add)), colnames(plain))), drop=FALSE])
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
    ### Schritt 6.1: das ".vs." durch minus ersetzen
           plain <- replaceVS(plain, out=out, groups=groups)
           colnames(plain) <- car::recode(colnames(plain), "'group'='label2'")
           if ( length(which(is.na(plain[,"comparison"]))) > 0) {plain[,"comparison"] <- car::recode(plain[,"comparison"], "NA='none'")}
    ### Schritt 7: post processing labelspalten 1 (sprechend) und 2 (pseudomathematisch) anpassen: 'p:\Methoden\02_R_Pakete\eatRep\Outputstruktur\Umfrage.doc'
           plain <- data.frame ( label1 = createLabel1(plain, out=out), plain, stringsAsFactors=FALSE)
browser()
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
           focCols<- unique(data.frame ( focCols, string[,c(attr(out, "allNam")[["trend"]], "depVar")],stringsAsFactors = FALSE))
           focus  <- merge(focCols,plainn, by = colnames(focCols), all=FALSE)
           if(nrow(focus) != length(unique(focus[,"parameter"]))) {focus  <- chooseLine(obj = focus, prepost=prepost, whichOne = 1, char = "_", groupName = "group")}
           if( prepost[[1]][[2]][1] == "wholeGroup") {
               refCol <- unique(data.frame ( group = "wholeGroup", string[,c(attr(out, "allNam")[["trend"]], "depVar")], stringsAsFactors = FALSE))
               refere <- merge(refCol,plainn, by = colnames(refCol), all=FALSE)
           }  else  {
               refCol <- data.frame ( lapply(plainn, FUN = function (col) { any(prepost[[1]][[2]] %in% col) })[attr(out, "allNam")[["group"]]], drop=FALSE, stringsAsFactors = FALSE)
               refCol <- refCol[,which(sapply(refCol, FUN = eval)), drop=FALSE]
               for ( i in 1:ncol(refCol)) {refCol[,i] <- prepost[[1]][[2]][i]}
               refCol <- unique(data.frame ( refCol, string[,c(attr(out, "allNam")[["trend"]], "depVar")], stringsAsFactors = FALSE))
               refere <- merge(refCol,plainn, by = colnames(refCol), all=FALSE)
               if ( length(prepost[[1]][[2]]) == 1) {refere  <- refere[setdiff(1:nrow(refere), grep("_", refere[,"group"])),]}
           }
           if(nrow(refere) != length(unique(refere[,"parameter"]))) {refere <- chooseLine(obj = refere, prepost=prepost, whichOne = 2, char = "_", groupName = "group")}
           return(rbind(focus, refere))}
           
chooseLine <- function(obj, prepost, whichOne, char, groupName) {
           anzahl <- stringr::str_count(obj[,groupName], char)
           needed <- length(prepost[[1]][[whichOne]]) - 1
           obj    <- obj[which(anzahl == needed),]
           return(obj)}

### Funktion waehlt die Fokus- und Referenzgruppe fuer trend differences aus
idFocRefGrpTREND <- function(string, out, plainn, add) {
           tomerge<- unique(string[,c(attr(out, "allNam")[["trend"]], "group", "depVar", names(add))])
           tomerge<- rbind(tomerge, tomerge)
           tomerge[,1] <- strsplit(tomerge[1,1], ".vs.")[[1]]
           tomerge<- merge(tomerge,plainn, by = colnames(tomerge), all=FALSE)
           stopifnot(nrow(tomerge)==2*length(unique(tomerge[,"parameter"])))
           return(tomerge)}

### Funktion waehlt die Fokus- und Referenzgruppe fuer group differences aus
idFocRefGrpGROUPDIFF <- function(string, out, plainn, add) {
           kontras<- unlist(lapply(unique(string[,attr(out, "allNam")[["group"]], drop=FALSE]), FUN = function (col) {grep(".vs.", col)}))
           stopifnot(length(kontras) == 1)
           kontr  <- strsplit(unique(string[,names(kontras)]), split = ".vs.")
           tomerge<- unique(string[,c(attr(out, "allNam")[["trend"]], "modus", "depVar", attr(out, "allNam")[["group"]], names(add))])
           tomerge<- rbind(tomerge, tomerge)
           tomerge[,names(kontras)] <- kontr[[1]]
           tomerge<- merge(tomerge,plainn, by = colnames(tomerge), all=FALSE)
           stopifnot(nrow(tomerge) == 2*length(unique(tomerge[,"parameter"])))
           return(tomerge)}

idFocRefGrpCROSSDIFF_OF_GROUPDIFF <- function (string, out, plainn, add) {
    ### erstmal den Gruppenvergleich finden (Referenz!)
           spl1   <- strsplit(unique(string[,"group"]), split = ".VS.")
           tomerge<- unique(string[,c(attr(out, "allNam")[["trend"]], "depVar", names(add))])
           tomerge[,"group"] <- spl1[[1]][2]
           refere <- merge(tomerge, subset(plainn, comparison == "groupDiff"), by = colnames(tomerge), all=FALSE)
           stopifnot(nrow(refere)==length(unique(refere[,"parameter"])))
    ### jetzt den crossdiff finden
           focus  <- unique(string[,c(attr(out, "allNam")[["trend"]], "depVar", names(add))])
           focus[,"group"] <- spl1[[1]][1]
           focus  <- merge(focus, subset(plainn, comparison == "groupDiff"), by = colnames(focus), all=FALSE)
           stopifnot(nrow(focus)==length(unique(focus[,"parameter"])))
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

setIDs <- function(comps, lev, plainc,oldNew,plainn, co1, indices, out, add) {
           if(nrow(plainc)>0 && nrow(subset(comps, level == lev))>0) {
               compsN<- do.call("rbind", apply(subset(comps, level == lev), MARGIN = 1, FUN = function (co) {
                        plc3 <- subset(plainc, comparison == co[[co1]])
                        plc3i<- do.call("rbind", by(data = plc3, INDICES = plc3[,indices], FUN = function (sg3) {
                                fun <- paste0("idFocRefGrp",eatTools::recodeLookup(toupper(unique(sg3[,"comparison"])), oldNew))
                                foc <- eval(parse(text=paste0(fun, "(string = sg3, out=out, plainn=plainn, add=add)")))
                                stopifnot(length(unique(foc[,"id"])) ==2)
                                sg3 <- data.frame ( sg3, id = paste("comp", sg3[1,"row"], sep="_"), unit_1 = sort(unique(foc[,"id"]))[1], unit_2 = sort(unique(foc[,"id"]))[2], stringsAsFactors = FALSE)
                                return(sg3)}))
                        return(plc3i)}, simplify = FALSE))
           }  else  {
               compsN<- NULL
           }
           return(compsN)}
           
replaceVS <- function (plain, out, groups) {
           if(length(c(attr(out, "allNam")[["group"]], attr(out, "allNam")[["trend"]]))>0) {for ( i in c(attr(out, "allNam")[["group"]], attr(out, "allNam")[["trend"]])) {
              if(i %in% colnames(groups)) {groups[,i] <- car::recode(groups[,i], "NA='total'")}
              if(i %in% colnames(plain)) {
                 plain[,i] <- car::recode(plain[,i], "NA='total'")
                 if(length(grep(".vs.",plain[,i]))>0) {
                    frme <- eatTools::halveString(plain[,i], ".vs.")
                    plain[,i] <- apply(X=frme, MARGIN=1, FUN = function (zeile) {
                                 if(length(na.omit(zeile)) ==1) {
                                    return(na.omit(zeile))
                                 } else {
                                    if ("total" %in% zeile) {
                                        return (paste0(zeile[-match("total", zeile)], " - total" ))
                                    } else {
                                        return (paste0(zeile[2], " - ", zeile[1]))
                                    }
                                 }})
                 }
              }
           }}
           return(plain)}
           
createLabel1 <- function(plain, out) {
           label1 <- apply(X=plain, MARGIN = 1, FUN = function (zeile){
                     if(!is.null(attr(out, "allNam")[["trend"]]) && grepl(" - ", zeile[[attr(out, "allNam")[["trend"]]]])) {
                        tr <- paste0("trend (",zeile[[attr(out, "allNam")[["trend"]]]],") fuer ")
                     } else {
                        tr <- ""
                     }
                     if(grepl("crossdiff", zeile[["comparison"]],ignore.case=TRUE)) {
                        cd <- lapply(attr(out, "allNam")[["group"]], FUN = function (v) {
                                spl <- unlist(strsplit(zeile[[v]], " - "))
                                if(length(spl)==1) {spl <- c(spl,spl)}
                                return(spl)})
                        cd <- paste(cd[[1]], cd[[2]], sep="_", collapse=" - ")
                        cd <- paste0("crossDiff (",cd,") ")
                     } else {
                        cd <- ""
                     }
                     if(grepl("groupdiff", zeile[["comparison"]],ignore.case=TRUE)) {
                        sig<- unlist(zeile[attr(out, "allNam")[["group"]]])
                        col<- setdiff(grep(" - ", sig), grep("total", sig))
                        stopifnot(length(col)==1)
                        gd <- paste0("bzgl. groupDiff (",sig[[col]],")")
                     } else {
                        gd <- ""
                     }
                     tot <- eatTools::crop(eatTools::crop(eatTools::crop(eatTools::crop(paste0(tr, cd, gd), "bzgl\\."), " "), "fuer"), " ")
                     return(tot)})
           return(label1)}
                     



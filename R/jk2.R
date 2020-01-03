generate.replicates <- function ( dat, ID, wgt = NULL, PSU, repInd, type )   {
                       if(type %in% c("JK2", "BRR")) { stopifnot(length(PSU) == 1 & length(repInd) == 1 ) }
                       if(type  == "JK1" ) { if(!is.null(repInd))  {
                          cat("'repInd' is ignored for 'type = JK1'.\n")
                          repInd <- NULL
                       }  }
                       allVars     <- list(ID = ID, wgt = wgt, PSU = PSU, repInd = repInd)
                       all.Names   <- lapply(allVars, FUN=function(ii) {existsBackgroundVariables(dat = dat, variable=ii)})
                       dat.i       <- dat[,unlist(all.Names)]
                       if(type %in% c("JK2", "BRR")) { if( !all( names(table(dat.i[,all.Names[["repInd"]]])) == c(0,1)) ) {stop("Only 0 and 1 are allowed for repInd variable.\n")} }
                       zonen       <- names(table(as.character(dat.i[,all.Names[["PSU"]]]) ) )
                       cat(paste("Create ",length(zonen)," replicate weights according to ",type," procedure.\n",sep=""))
                       if ( nrow(dat)>2500 & length(zonen) > 50 ) {             ### progress bar wird nur angezeigt, wenn es auch nennenswerte Zeit dauert
                            pb     <- progress_bar$new( format = "  replicates [:bar] :percent in :elapsed", incomplete = " ", total = length(zonen), clear = FALSE, width= 60, show_after = 0.01)
                       }
                       missings    <- sapply(dat.i, FUN = function (ii) {length(which(is.na(ii)))})
                       if(!all(missings == 0)) {
                           mis.vars <- paste(names(missings)[which(missings != 0)], collapse = ", ")
                           stop(paste("Found missing value(s) in variable(s) ", mis.vars,".\n",sep=""))
                       }                                                        ### untere Zeile: in JK-2 werden die gewichte nur in der entsprechenden PSU _fuer ein Replicate_ geaendert
                       reps <- data.frame ( lapply(zonen , FUN = function(ii) { ### in BRR werden die Gewichte in allen Zonen _fuer ein Replicate_ geaendert!
                               if ( nrow(dat)>2500 & length(zonen) > 50 ) { pb$tick() }
                               rep.ii <- dat.i[,all.Names[["wgt"]]]             ### in JK-1 werden die Gewichte nur in der entsprechenden PSU _fuer alle Replicates_ geaendert
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

### Wrapper: ruft "eatRep()" mit selektiven Argumenten auf
jk2.mean <- function(datL, ID, wgt = NULL, type = c("JK2", "JK1", "BRR"), PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL,
            group.splits = length(groups), group.differences.by = NULL, cross.differences = FALSE, crossDiffSE = c("wec", "pisa","old"), nBoot = 100,
            group.delimiter = "_", trend = NULL, linkErr = NULL, dependent, na.rm = FALSE, doCheck = TRUE, engine = c("survey", "BIFIEsurvey") ) {
            cdse<- match.arg(arg = crossDiffSE, choices = c("wec", "pisa","old"))
            type<-  match.arg(arg = toupper(type), choices = c("JK2", "JK1", "BRR"))
            if(!"data.frame" %in% class(datL) || "tbl" %in% class(datL) ) { cat(paste0("Convert 'datL' of class '",paste(class(datL), collapse="', '"),"'to a data.frame.\n")); datL <- data.frame ( datL, stringsAsFactors = FALSE)}
            if ( is.null ( attr(datL, "modus"))) {
                  modus <- identifyMode ( name = "mean", type = type, PSU = PSU, repWgt=repWgt )
            }  else  {
                  modus <- attr(datL, "modus")
            }
            ret <- eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "mean",
                   nest = nest, imp = imp, groups = groups, group.splits = group.splits, group.differences.by = group.differences.by,
                   cross.differences = cross.differences, trend = trend, linkErr = linkErr, dependent = dependent, group.delimiter=group.delimiter, na.rm=na.rm, doCheck=doCheck, modus = modus, engine=engine)
            if ( isFALSE(ret[["allNam"]][["cross.differences"]])) {return(ret)}
            if ( (is.list(cross.differences) || cross.differences == TRUE) && cdse == "old" ) {
                 ret[["SE_correction"]] <- list(NULL)
                 class(ret[["SE_correction"]]) <- c("old", "list")
                 return(ret)
            }
            if ( is.list(cross.differences) || isTRUE(cross.differences) ) {
                  toAppl<- eatRep:::superSplitter(group = ret[["allNam"]][["group"]], group.splits = group.splits, group.differences.by = ret[["allNam"]][["group.differences.by"]], group.delimiter = group.delimiter , dependent=ret[["allNam"]][["dependent"]] )
                  stopifnot(length ( toAppl ) > 1)
                  vgl <- do.call("rbind", lapply(1:length(toAppl), FUN = function ( y ) {
                         gdb <- attr(toAppl[[y]], "group.differences.by")
                         if ( is.null(gdb)) {gdb <- NA}
                         res <-  data.frame ( analysis.number = y, hierarchy.level = length(toAppl[[y]]), groups.divided.by = paste(toAppl[[y]], collapse=" + "), group.differences.by = gdb)
                         return(res)}))
                  ana <- lapply ( combn(x=vgl[,"analysis.number"], m=2, simplify=FALSE), FUN = function ( a ) {
                         t1 <- abs(diff(vgl[ match(a, vgl[,"analysis.number"]),"hierarchy.level"])) == 1
                         t2 <- TRUE                                             ### Vergleiche sollen nur angestellt werden, wenn Differenz der Hierarchielevel = 1 und
                         if ( vgl[min(a),"hierarchy.level"] != 0 ) {            ### alle bis auf eine Gruppierungsvariable der hoeheren Ebene in der unteren enthalten sind, und
                              nam<- lapply(sort(a), FUN = function ( z ) {      ### wenn Vergleiche Teilmenge dessen ist, was in 'cross.differences' angegeben ist
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
                      cat(paste0("Warning: Computation of cross level differences using '",cdse,"' method is only possible for differences according to adjacent levels. Non-adjacent levels will be ignored.\n"))
                  }  else {
                      cat(paste0("Compute cross level differences using '",cdse,"' method.\n"))
                  }
                  spl <- lapply(ana, FUN = function ( a ) {
                         vgl <- vgl[a,]
                         if (vgl[1,"hierarchy.level"] != 0) {                   ### wenn referenzdatensatz bspw. alle Maedchen sein soll,wird der jetzt nach den Gruppierungsvariablen gesplittet
                             nam <- unlist(strsplit(as.character(vgl[1,"groups.divided.by"]), split=" |\\+"))
                             nam <- nam[which(nchar(nam)>0)]
                             dat <- by(datL, INDICES = datL[,nam], FUN = function ( d ) {return(d)})
                             grp <- setdiff(unlist(strsplit(as.character(vgl[2,"groups.divided.by"]), split=" |\\+")) , nam)
                             grp <- grp[which(nchar(grp)>0)]
                         }  else  {                                             ### else = Referenz ist Gesamtdatensatz -> Liste mit nur einem Element
                             dat <- list(datL)
                             grp <- as.character(vgl[2,"groups.divided.by"])
                         }
                         stopifnot(length(grp)==1)
                         cld <- lapply(dat, FUN = function ( d ) {              ### 'cld' = cross level differences
                                if ( is.null(ret[["allNam"]][["wgt"]]) || !ret[["allNam"]][["wgt"]] %in% colnames(d)) {
                                     stopifnot(is.null(wgt))                    ### workaround: wenn keine gewichte spezifiziert wurden, werden sie auf 1 gesetzt und allNam$wgt ist nicht NULL
                                     cat(paste0("   '",cdse,"' method: Assume equally weighted cases.\n"))
                                     gew <- NULL                                ### deswegen sucht er im Datensatz eine gewichtungsvariable, die nur 1en enthaelt und gibt Fehlermeldung
                                } else {                                        ### Fuer diesen Fall wird fuer 'wgt' das Argument 'NULL' uebergeben
                                     gew <- ret[["allNam"]][["wgt"]]
                                }                                               ### untere Zeilen: bloeder hotfix, notwendig wegen anderem Hotfix, weil ret$allNam$nest nicht NULL ist, selbst wenn keine nests spezifiziert werden (damit 'by' funktioniert)
                                if (is.null(nest)) {ne <- NULL} else {ne <- ret[["allNam"]][["nest"]]}
                                if (is.null(imp))  {im <- NULL} else {im <- ret[["allNam"]][["imp"]]}
                                if ( vgl[1,"hierarchy.level"] != 0) {
                                     rg <- facToChar(d[1,nam, drop=FALSE])
                                     rg <- do.call("rbind", lapply(names(rg), FUN = function (r){data.frame ( groupName = r, groupValue = rg[1,r], stringsAsFactors = FALSE) }))
                                }  else  {
                                     rg <- NULL
                                }
                                if ( cdse == "wec" ) {
                                     b <- jk2.glm(datL=d, ID=ret[["allNam"]][["ID"]], wgt = gew, type = type, PSU = ret[["allNam"]][["PSU"]], repInd = ret[["allNam"]][["repInd"]],
                                                  repWgt = ret[["allNam"]][["repWgt"]], nest=ne, imp=im, trend = trend,
                                                  formula = as.formula(paste0(ret[["allNam"]][["dependent"]] , " ~ ", grp)), doCheck = doCheck, na.rm = na.rm, useWec = TRUE, engine = engine )
                                }  else  {                                      ### hier beginnt methode 'pisa'
                                     b <- eatRep(datL =d, ID=ret[["allNam"]][["ID"]], wgt = gew, type=type, PSU = ret[["allNam"]][["PSU"]], repInd = ret[["allNam"]][["repInd"]], toCall = "cov",nBoot=nBoot,
                                          nest=ne, imp=im, groups = grp, refGrp = rg, trend = trend, dependent = ret[["allNam"]][["dependent"]], na.rm=na.rm, doCheck=FALSE, engine=engine, modus=modus)
                                }
                                b[["vgl"]]      <- vgl
                                b[["focGrp"]]   <- grp                          ### Fokus- und Referenzgruppe in Ergebnisoutput eintragen,
                                b[["refGrp"]]   <- "all"                        ### zur Weiterverarbeitung in den Reportingfunktionen
                                if(!is.null(rg)) {b[["refGrp"]] <- rg}
                         return(b)})
                  return(cld)})
                  spl <- unlist(spl, recursive=FALSE)
                  class(spl) <- c( recode(cdse, "'wec'='wec_se_correction'; 'pisa'='pisa_se_correction'"), "list")
                  ret[["SE_correction"]] <- spl
            }
            return(ret)}


### Wrapper: ruft "eatRep()" mit selektiven Argumenten auf
jk2.table<- function(datL, ID, wgt = NULL, type = c("JK2", "JK1", "BRR"),
            PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL, refGrp = NULL, group.splits = length(groups), group.differences.by = NULL, cross.differences = FALSE, crossDiffSE = c("wec", "pisa","old"),
            nBoot = 100, chiSquare = FALSE, correct = TRUE, group.delimiter = "_", trend = NULL, linkErr = NULL, dependent , separate.missing.indicator = FALSE,na.rm=FALSE, expected.values = NULL, doCheck = TRUE, forceTable = FALSE,
            engine = c("survey", "BIFIEsurvey") ) {                             ### untere Zeile: wrapper! hier wird jk2.table ueber jk2.mean aufgerufen; fuer eine kategorielle Variable sind die Haeufigkeiten die Mittelwerte der Indikatoren der Faktorstufen; untere Zeilen, Achtung!! hier muessen immer zwei '&'-Zeichen gesetzt werden!!
            modus <- identifyMode ( name = "table", type = type, PSU = PSU, repWgt=repWgt )
           if(!"data.frame" %in% class(datL) || "tbl" %in% class(datL) ) { cat(paste0("Convert 'datL' of class '",paste(class(datL), collapse="', '"),"'to a data.frame.\n")); datL <- data.frame ( datL, stringsAsFactors = FALSE)}
            chk1  <- eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "table",
                     nest = nest, imp = imp, groups = groups, group.splits = group.splits, group.differences.by = group.differences.by, cross.differences=cross.differences, correct = correct,
                     trend = trend, linkErr = linkErr, dependent = dependent, group.delimiter=group.delimiter, separate.missing.indicator=separate.missing.indicator,
                     expected.values=expected.values, na.rm=na.rm, doCheck=FALSE, onlyCheck= TRUE, modus = modus, engine=engine)
            if ( length(unique(datL[,chk1[["dependent"]]])) == 2 && all(sort(unique(datL[,chk1[["dependent"]]])) == 0:1) == TRUE && forceTable == FALSE ) {
                 attr(datL, "modus") <- modus
                 ret <- jk2.mean ( datL = datL, ID=chk1[["ID"]], wgt=chk1[["wgt"]], type = type, PSU = chk1[["PSU"]], repInd = chk1[["repInd"]], repWgt = repWgt,
                        nest = chk1[["nest"]], imp = chk1[["imp"]], groups = groups, group.splits = group.splits, group.differences.by=group.differences.by, cross.differences=cross.differences,
                        crossDiffSE = crossDiffSE, nBoot = nBoot, group.delimiter =group.delimiter, trend = chk1[["trend"]], linkErr = chk1[["linkErr"]], dependent = chk1[["dependent"]], na.rm=na.rm,doCheck = doCheck, engine=engine)
                 return(ret)
            }  else  {
                 if ( !is.null(group.differences.by) && isFALSE(chiSquare)) {
                    chk <- eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "table",
                           nest = nest, imp = imp, groups = groups, group.splits = group.splits, group.differences.by = group.differences.by, cross.differences=cross.differences, correct = correct,
                           trend = trend, linkErr = linkErr, dependent = dependent, group.delimiter=group.delimiter, separate.missing.indicator=separate.missing.indicator,
                           expected.values=expected.values, na.rm=na.rm, doCheck=doCheck, onlyCheck= TRUE, modus = modus, engine=engine)
    ### missing handling muss vorneweg geschehen
                    isNa<- which ( is.na ( datL[, chk[["dependent"]] ] ))
                    if ( length ( isNa ) > 0 ) {
                         cat(paste0("Warning: Found ",length(isNa)," missing values in dependent variable '",chk[["dependent"]],"'.\n"))
                         if ( isTRUE(separate.missing.indicator) ) {
                              stopifnot ( length( intersect ( "missing" , names(table(datL[, chk[["dependent"]] ])) )) == 0 )
                              if(class(datL[, chk[["dependent"]] ]) == "factor") {# Hotfix: fuer Faktorvariablen funktioniert das einfache subsetting 
                                  levOld <- levels(datL[, chk[["dependent"]] ]) ### dat[which(is.na(dat[,"var"])) ,"var"] <- "missing" nicht
                                  datL[, chk[["dependent"]] ] <- as.character(datL[, chk[["dependent"]] ])
                                  datL[isNa, chk[["dependent"]] ] <- "missing"
                                  datL[, chk[["dependent"]] ] <- factor(datL[, chk[["dependent"]] ], levels=c(levOld, "missing"))
                              }  else  {
                                  datL[isNa, chk[["dependent"]] ] <- "missing"
                              }
                         }  else  {
                              if ( na.rm == FALSE ) { stop("If no separate missing indicator is used ('separate.missing.indicator == FALSE'), 'na.rm' must be TRUE if missing values occur.\n")}
                              datL <- datL[-isNa,]
                         }
                    }                                                           ### Ende des missing handlings
                    frml<- as.formula ( paste("~ ",chk[["dependent"]]," - 1",sep="") )
                    datL[, chk[["dependent"]] ] <- as.character( datL[, chk[["dependent"]] ] )
                    matr<- data.frame ( model.matrix ( frml, data = datL) )     ### untere zeilen: Wrapper liefert objekt von jk2.mean zurueck,
                    datL<- data.frame ( datL,  matr)                            ### der Output muss deshalb in das Format von 'table' umgeformt werden, das macht die Funktion 'clearTab' am Ende von 'eatRep'
                    ret <- lapply ( colnames(matr), FUN = function ( dpd ) {
                           attr(datL, "modus") <- modus
                           attr(datL,"depOri") <- chk[["dependent"]]
                           res <- jk2.mean ( datL = datL, ID=chk[["ID"]], wgt=chk[["wgt"]], type = type, PSU = chk[["PSU"]], repInd = chk[["repInd"]], repWgt = repWgt,
                                             nest = chk[["nest"]], imp = chk[["imp"]], groups = groups, group.splits = group.splits, group.differences.by=group.differences.by, cross.differences=cross.differences,
                                             crossDiffSE = crossDiffSE, nBoot = nBoot, group.delimiter =group.delimiter, trend = chk[["trend"]], linkErr = chk[["linkErr"]], dependent = dpd, na.rm=na.rm,doCheck = doCheck, engine=engine)
                           return(res) } )
    ### Output zusammenfassen: kein Trend ... dann ist "resT" eine Liste mit nur einem Element, das ein data.frame ist
                    if ( length(ret[[1]][["resT"]]) == 1) {
                         lst <- do.call("rbind", lapply(ret, FUN = function ( x ) { x[["resT"]][["noTrend"]]}))
                         ret <- ret[[1]]
                         ret[["resT"]][["noTrend"]] <- lst
                    }  else  {
    ### Output zusammenfassen: Trend
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
                           expected.values=expected.values, na.rm=na.rm, doCheck=doCheck, modus=modus, engine=engine)
                 }
                 return(ret)}}


### Wrapper: ruft "eatRep()" mit selektiven Argumenten auf
jk2.quantile<- function(datL, ID, wgt = NULL, type = c("JK2", "JK1", "BRR"),
            PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL, group.splits = length(groups),
            cross.differences = FALSE, group.delimiter = "_", trend = NULL, linkErr = NULL, dependent, probs = seq(0, 1, 0.25),  na.rm = FALSE,
            nBoot = NULL, bootMethod = c("wSampling","wQuantiles") , doCheck = TRUE, engine = c("survey", "BIFIEsurvey"))  {
            modus      <- identifyMode ( name = "quantile", type = type, PSU = PSU, repWgt=repWgt )
            bootMethod <- match.arg ( bootMethod )
            eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "quantile",
                   nest = nest, imp = imp, groups = groups, group.splits = group.splits, cross.differences=cross.differences, trend = trend, linkErr = linkErr, dependent = dependent,
                   group.delimiter=group.delimiter, probs=probs, na.rm=na.rm, nBoot=nBoot, bootMethod=bootMethod, doCheck=doCheck, modus=modus, engine=engine)}


### Wrapper: ruft "eatRep()" mit selektiven Argumenten auf
jk2.glm  <- function(datL, ID, wgt = NULL, type = c("JK2", "JK1", "BRR"),
            PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL, group.splits = length(groups), group.delimiter = "_", cross.differences = FALSE,
            trend = NULL, linkErr = NULL, formula, family=gaussian, forceSingularityTreatment = FALSE, glmTransformation = c("none", "sdY"),
            doCheck = TRUE, na.rm = FALSE, poolMethod = c("mice", "scalar") , useWec = FALSE, engine = c("survey", "BIFIEsurvey")) {
            modus  <- identifyMode ( name = "glm", type = type, PSU = PSU, repWgt=repWgt )
            poolMethod <- match.arg(poolMethod)
            eatRep(datL =datL, ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "glm",
                   nest = nest, imp = imp, groups = groups, group.splits = group.splits, cross.differences = cross.differences, trend = trend, linkErr = linkErr,
                   formula=formula, family=family, forceSingularityTreatment=forceSingularityTreatment, glmTransformation = glmTransformation,
                   group.delimiter=group.delimiter, na.rm=na.rm, doCheck=doCheck, modus=modus, poolMethod=poolMethod, useWec=useWec, engine=engine)}


### Funktion ist nicht user-level, sondern wird von jk2.mean, jk2.table, jk2.quantile, jk2.glm mit entsprechenden Argumenten aufgerufen
eatRep <- function (datL, ID, wgt = NULL, type = c("JK2", "JK1", "BRR"), PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL,
          toCall = c("mean", "table", "quantile", "glm", "cov"), groups = NULL, refGrp = NULL, group.splits = length(groups), group.differences.by = NULL,
          cross.differences = FALSE, group.delimiter = "_", trend = NULL, linkErr = NULL, dependent, na.rm = FALSE, forcePooling = TRUE, boundary = 3, doCheck = TRUE,
          separate.missing.indicator = FALSE, expected.values = NULL, probs = NULL, nBoot = NULL, bootMethod = NULL, formula=NULL, family=NULL,
          forceSingularityTreatment = FALSE, glmTransformation = c("none", "sdY"), correct, onlyCheck = FALSE, modus, poolMethod = "mice", useWec = FALSE, engine) {
          if(!"data.frame" %in% class(datL) || "tbl" %in% class(datL) ) { cat(paste0("Convert 'datL' of class '",paste(class(datL), collapse="', '"),"'to a data.frame.\n")); datL <- data.frame ( datL, stringsAsFactors = FALSE)}
          if ( useWec == TRUE ) { forceSingularityTreatment <- TRUE; poolMethod <- "scalar"}
          i     <- 0                                                            ### wenn wec genutzt werden soll, geht das nur ueber den Hotfix, der auch fuer Singularitaet verwendet wird
          fc    <- NULL                                                         ### initialisieren
          while ( !is.null(sys.call(i))) { fc <- c(fc, crop(unlist(strsplit(deparse(sys.call(i))[1], split = "\\("))[1])); i <- i-1  }
          fc   <- fc[max(which(fc %in% c("jk2.mean", "jk2.table", "jk2.glm", "jk2.quantile")))]
          toCall<- match.arg(toCall)                                            ### 'oberste' Funktion suchen, die eatRep gecallt hat; zweiter Teil des Aufrufs ist dazu da, dass nicht "by" drinsteht, wenn "jk2.mean" innerhalb einer anderen "by"-Funktion aufgerufen wird
          type  <- match.arg(arg = toupper(type), choices = c("JK2", "JK1", "BRR"))
          engine<- match.arg(arg = engine, choices = c("survey", "BIFIEsurvey"))
          glmTransformation <- match.arg(glmTransformation)
          if(forceSingularityTreatment == FALSE & glmTransformation != "none") {
             cat("'forceSingularityTreatment' was set to 'FALSE'. Please note that 'glmTransformation' is only possible if 'forceSingularityTreatment' is 'TRUE'.\n"); flush.console()
          }
          if(toCall == "glm") {                                                 ### fuer glm muessen abhaengge und unabhaengige Variablen aus Formel extrahiert werden
             dependent  <- as.character(formula)[2]
             independent<- all.vars(formula[-2])
           }  else {
             independent <- NULL
          }
          if(is.null(groups))  {groups <- "wholeGroup"; datL[,"wholeGroup"] <- "wholePop"}
          allVar<- list(ID = ID, wgt = wgt, PSU = PSU, repInd = repInd, repWgt = repWgt, nest=nest, imp=imp, group = groups, trend=trend, linkErr = linkErr, group.differences.by=group.differences.by, dependent = dependent, independent=independent)
          allNam<- lapply(allVar, FUN=function(ii) {existsBackgroundVariables(dat = datL, variable=ii)})
          if(forceSingularityTreatment == TRUE && !is.null(allNam[["PSU"]]) ) { poolMethod <- "scalar"}
    ### check: Gruppierungsvariablen duerfen nicht numerisch sein
          if(!is.null(allNam[["group"]])) {
             chk <- lapply(allNam[["group"]], FUN = function ( v ) { if ( !class(datL[,v]) %in% c("factor", "character", "logical", "integer")) {stop(paste0("Grouping variable '",v,"' must be of class 'factor', 'character', 'logical', or 'integer'.\n"))} })
             if ( is.list(cross.differences) || isTRUE(cross.differences)) {    ### levels der Gruppen duerfen keine "." oder "_" enthalten, falls cross differences berechnet werden sollen
                 for ( gg in allNam[["group"]] ) {
                       if ( length(grep("\\.|_", datL[,gg])) > 0) {
                           cat(paste0( "Levels of grouping variable '",gg, "' contain '.' and/or '_' which is not allowed. '.' and '_' will be deleted.\n"))
                           if ( class ( datL[,gg] ) == "factor") {
                               levNew <- gsub("\\.|_", "", levels(datL[,gg]))
                               datL[,gg] <- factor(gsub("\\.|_", "", datL[,gg]), levels = levNew)
                           }  else  {
                               datL[,gg] <- gsub("\\.|_", "", datL[,gg])
                           }
                       }
                 }
             }
          }
    ### check: Manche Variablennamen sind in der Ergebnisstruktur fest vergeben und duerfen daher nutzerseitig nicht als namen von Gruppenvariablen vergeben werden
          na    <- c("isClear", "N_weightedValid", "N_weighted",  "wgtOne")
          naGr  <- c("wholePop", "group", "depVar", "modus", "parameter", "coefficient", "value", "linkErr", "comparison", "sum", "trendvariable")
          naInd <- c("(Intercept)", "Ncases", "Nvalid", "R2",  "R2nagel", "linkErr")
          naGr1 <- which ( allNam[["group"]] %in% naGr )                        ### hier kuenftig besser: "verbotene" Variablennamen sollen automatisch umbenannt werden!
          if(length(naGr1)>0)  {stop(paste0("Following name(s) of grouping variables in data set are forbidden due to danger of confusion with result structure:\n     '", paste(allNam[["group"]][naGr1], collapse="', '"), "'\n  Please rename these variable(s) in the data set.\n"))}
          naInd1<- which ( allNam[["independent"]] %in% naInd )
          if(length(naInd1)>0)  {stop(paste0("Following name(s) of independent variables in data set are forbidden due to danger of confusion with result structure:\n     '", paste(allNam[["independent"]][naInd1], collapse="', '"), "'\n  Please rename these variable(s) in the data set.\n"))}
          na2   <- which ( unlist(allNam) %in% na )
          if(length(na2)>0)  {stop(paste0("Following variable name(s) in data set are forbidden due to danger of confusion with result structure:\n     '", paste(unlist(allNam)[na2], collapse="', '"), "'\n  Please rename these variable(s) in the data set.\n"))}
    ### fuer weighted effect coding muss UV ein Faktor sein, und es darf auch nur eine UV geben
          if (isTRUE(useWec) ) {
              if ( length(independent) != 1 ) {stop("Only one independent (grouping) variable is allowed for weighted effect coding.\n")}
              if ( !class(datL[,independent]) %in% c("factor", "character", "logical", "integer")) {stop(paste0("For weighted effect coding, independent (grouping) variable '",independent,"' must be of class 'factor', 'character', 'logical', or 'integer'.\n"))}
          }
    ### Methode "mice" geht nur fuer nicht-genestete Imputationen
          if (poolMethod == "mice" && !is.null(allNam[["nest"]]) && toCall == "glm" ) {
              cat("Method 'mice' is not available for nested imputation. Switch to method 'scalar'.\n")
              poolMethod <- "scalar"
          }
    ### Achtung: wenn Funktion nur zum checken genutzt wird, endet sie hier
          if ( isTRUE(onlyCheck) ) {
              ret <- allNam
          }  else  {
    ### wie in 'defineModel': Funktion ruft sich rekursiv selber auf, wenn Trend bestimmt werden soll
              if(!is.null(allNam[["trend"]])) {                                 ### Achtung: hier wenn Trendberechnung geschehen soll
    ### check: nur zwei Gruppen in Trendvariablen?
                  lev <- sort ( unique(datL[,allNam[["trend"]]]))
                  if(length(lev) != 2) {stop(paste(length(lev), " levels ('",paste(lev, collapse="', '"),"') found for the 'trend' variable '",allNam[["trend"]],"'. 2 levels are allowed.\n",sep=""))}
    ### check: alle Kombinationen von faktoren in beiden datensaetzen?
                  if (!is.null(allNam[["group"]])) {
                       foo <- lapply(allNam[["group"]], FUN = function ( gr ) {
                              ch <- by(data = datL, INDICES = datL[,allNam[["trend"]]], FUN = function ( subdat ) { table(subdat[,gr]) }, simplify = FALSE )
                              if ( !all ( names(ch[[1]]) == names(ch[[2]]))) { stop(paste("Error in grouping variable '",gr,"': Levels do not match. Levels in trend group '",names(ch)[1],"': \n    ", paste(names(ch[[1]]), collapse = ", "),"\nLevels in trend group '",names(ch)[2],"': \n    ", paste(names(ch[[2]]), collapse = ", "),sep="")) } } )
                  }
    ### repliziere Analyse zweimal, d.h. fuer beide Trendgruppen ... erstmal single core ... Funktion ruft sich 2x selber auf
                  resT<- by ( data = datL, INDICES = datL[,allNam[["trend"]]], FUN = function ( subdat ) {
                         cat(paste("\nTrend group: '",subdat[1,allNam[["trend"]] ], "'\n",sep=""))
                         do    <- paste ( "foo <- eatRep ( ", paste(names(formals(eatRep)), recode(names(formals(eatRep)), "'trend'='NULL'; 'datL'='subdat'; 'group.differences.by'='group.differences.by'"), sep =" = ", collapse = ", "), ")",sep="")
                         eval(parse(text=do))
                         foo[["resT"]][["noTrend"]][,allNam[["trend"]]] <- subdat[1,allNam[["trend"]] ]
                         out1  <- foo[["resT"]][["noTrend"]]
                         out2  <- foo[["allNam"]]
                         return(list(out1=out1, out2=out2))}, simplify = FALSE)
    ### Trends wurden berechnet: 'resT' aufbereiten, Linkingfehler ranmatchen ... wenn es keinen linkingfehler gibt, wird er auf 0 gesetzt
                  allNam <- c(allNam, resT[[1]][["out2"]])
                  allNam <- allNam[unique(names(allNam))]
                  le     <- createLinkingError ( allNam = allNam, resT = resT, datL = datL, fc=fc, toCall = toCall)
                  out3   <- lapply(resT, FUN = function ( k ) { k[["out1"]]})
                  ret    <- list(resT = out3, allNam = allNam, toCall = toCall, family=family, le=le)
                  return(ret)
              }  else {
    ### obere Zeile: Ende der inneren Schleife (= Ende des Selbstaufrufs wegen trend) ... das untere wird nun fuer jeden Aufruf abgearbeitet
                  if( length( setdiff ( allNam[["group.differences.by"]],allNam[["group"]])) != 0) {stop("Variable in 'group.differences.by' must be included in 'groups'.\n")}
                  if(toCall == "glm") {.GlobalEnv$glm.family <- family}         ### Hotfix!
    ### Anzahl der Analysen aufgrund mehrerer Hierarchieebenen definieren ueber den 'super splitter' und Analysen einzeln (ueber 'lapply') starten
                  toAppl<- superSplitter(group = allNam[["group"]], group.splits = group.splits, group.differences.by = allNam[["group.differences.by"]], group.delimiter = group.delimiter , dependent=allNam[["dependent"]] )
                  cat(paste(length(toAppl)," analyse(s) overall according to: 'group.splits = ",paste(group.splits, collapse = " ") ,"'.", sep=""))
    ### bei mehr als einer Analyse: genauere Ausgabe auf Konsole printen und consistency checks
                  if ( length ( toAppl ) > 1) {
                       ret <- do.call("rbind", lapply(1:length(toAppl), FUN = function ( y ) {
                              gdb <- attr(toAppl[[y]], "group.differences.by")
                              if ( is.null(gdb)) {gdb <- NA}
                              res <-  data.frame ( analysis.number = y, hierarchy.level = length(toAppl[[y]]), groups.divided.by = paste(toAppl[[y]], collapse=" + "), group.differences.by = gdb)
                              return(res)}))
                       cat("\n \n"); print(ret)
                  }
    ### cross differences defaulten wenn vom user nicht definiert
                  if ( !is.list(cross.differences)) {
                       if ( !is.logical (cross.differences)) {
                            cat("Argument 'cross.differences' must be either logical or a list. Set 'cross.differences' to FALSE.\n")
                            cross.differences <- FALSE
                       }
                       if ( isTRUE(cross.differences)) {
                            if ( "wholeGroup" %in% allNam[["group"]] ) {
                                  cat("No groups defined. Set 'cross.differences' to FALSE.\n")
                                  cross.differences <- FALSE
                            }  else  {
                                  if ( length(group.splits)==1) {
                                        cat("\nArgument 'group.splits' was set to 1. No cross-level differences can be computed. Set 'cross.differences' to FALSE.\n")
                                        cross.differences <- FALSE
                                  }  else  {
                                        cross.differences <- combn(x=group.splits,m=2, simplify = FALSE)
                                  }
                            }
                       }
                  }  else  {                                                    ### checks muessen nur gemacht werden, wenn die Liste nicht automatisch erzeugt wurde
                       if(!all(unlist(lapply(cross.differences, length)) == 2)) {stop("Each element in 'cross.differences' must be a vector of length 2.\n")}
                       if(!all(unlist(lapply(cross.differences, class)) %in% c("integer","numeric"))) {stop("Each element in 'cross.differences' must be a numerical vector.\n")}
                       if(!all(unlist(lapply(cross.differences, FUN = function ( x ) { length(unique(x))})) ==2 ) ) {stop("Each element in 'cross.differences' must be a vector of 2 different numbers.\n")}
                       vals <- unique(unlist(cross.differences))
                       if(!all(vals %in% group.splits)) {stop("All numerical values in 'cross.differences' must be included in 'group.splits'.\n")}
                       cross.differences <- lapply( cross.differences, sort )
                       if ( length(cross.differences) != length(unique(cross.differences)) ) {
                            cat("Some comparisons in 'cross.differences' are identical. Duplicated comparisons will be removed.\n")
                            cross.differences <- unique(cross.differences)
                       }
                  }
                  allNam[["cross.differences"]] <- cross.differences
    ### Achtung: wenn keine Gruppen und/oder Nests und/oder Imputationen spezifiziert sind, erzeuge Variablen mit Werten gleich 1, damit by() funktioniert!
                  if( is.null(allNam[["imp"]]) )  { datL[,"imp"] <- 1; allNam[["imp"]] <- "imp" } else { stopifnot(length(allNam[["imp"]]) == 1 ); datL[,allNam[["imp"]]] <- as.character(datL[,allNam[["imp"]]])}
                  if( is.null(allNam[["wgt"]]) )  { datL[,"wgtOne"] <- 1; allNam[["wgt"]] <- "wgtOne" } else {
                      stopifnot(length(allNam[["wgt"]]) == 1 )
                      if ( !class(datL[,allNam[["wgt"]]]) %in% c("numeric", "integer") ) { stop ( paste("Error: 'wgt' variable '",allNam[["wgt"]],"' of class '",class(datL[,allNam[["wgt"]]]),"' has to be numeric.\n",sep="")) }
                      isMis <- which(is.na(datL[,allNam[["wgt"]]]))
                      isZero<- which ( datL[,allNam[["wgt"]]] == 0 )
                      if(length(isMis)>0) { stop (paste ( "Error: Found ",length(isMis)," missing values in the weight variable '",allNam[["wgt"]],"'.\n",sep="")) }
                      if(length(isZero)>0) { cat (paste ( "\nWarning: Found ",length(isZero)," zero weights in the weight variable '",allNam[["wgt"]],"'.\n",sep="")) }
                  }
                  if(!is.null(allNam[["nest"]]))  {
                      stopifnot(length(allNam[["nest"]]) == 1 )
                      datL[,allNam[["nest"]]] <- as.character(datL[,allNam[["nest"]]])
                      cat(paste("\nAssume nested structure with ", length(table(datL[,allNam[["nest"]]]))," nests and ",length(table(datL[,allNam[["imp"]]]))," imputations in each nest. This will result in ",length(table(datL[,allNam[["nest"]]]))," x ",length(table(datL[,allNam[["imp"]]]))," = ",length(table(datL[,allNam[["nest"]]]))*length(table(datL[,allNam[["imp"]]]))," imputation replicates.\n",sep=""))
                  }  else  { cat("\nAssume unnested structure with ",length(table(datL[,allNam[["imp"]]]))," imputations.\n",sep="")}
                  datL[,"isClear"] <- TRUE
                  if( is.null(allNam[["nest"]]) ) { datL[,"nest"]  <- 1; allNam[["nest"]]  <- "nest" }
                  if(!is.null(allNam[["group"]])) {                             ### untere Zeile: das, damit leere Gruppen nicht ueber by() mit geschleift werden, wie es passiert, wenn Gruppen als Faktoren definiert sind
                      for ( jj in allNam[["group"]] )  { datL[,jj] <- as.character(datL[,jj]) }
                  }
    ### check: wenn cross.differences gemacht werden sollen, dann muessen die faktor levels aller gruppierungsvariablen disjunkt sein (siehe Mail benjamin, 13.11.2019, 18.11 Uhr)
                  if(!is.null(allNam[["cross.differences"]])) {
                      if(length(allNam[["group"]])>1) {
                         lev <- unlist(lapply(allNam[["group"]], FUN = function ( v ) { unique(as.character(datL[,v]))}))
                         if (length(lev) != length(unique(lev))) {stop("Factor levels of grouping variables are not disjunct.\n")}
                      }
                  }
    ### check: abhaengige Var. numerisch?
                  if(toCall %in% c("mean", "quantile", "glm")) {
                     if(!class(datL[,allNam[["dependent"]]]) %in% c("integer", "numeric")) {
                         stop(paste0("Dependent variable '",allNam[["dependent"]],"' has to be of class 'integer' or 'numeric'.\n"))
                     }
                  }
    ### wenn Replicates bereits uebergeben, muss PSU und repInd NULL sein
                  if(!is.null(repWgt) ) {
                     if ( !is.null(allNam[["PSU"]]) | !is.null(allNam[["repInd"]]) ) {
                         cat("Warning: Arguments 'PSU' and 'repInd' are expected to be NULL if replicate weights are already defined (via 'repWgt').\n    'PSU' and 'repInd' will be ignored.\n")
                     }
                  }
    ### replicates erzeugen (nur einmal fuer alle Analysen, und nur wenn 'repWgt' NULL ist)
                  if(!is.null(allNam[["repWgt"]]))  {
                     repA <- data.frame ( datL[,allNam[["ID"]], drop=FALSE], datL[,allNam[["repWgt"]] ])
                     repA <- repA[!duplicated(repA[,allNam[["ID"]]]),]
                  }  else  {
                     if(!is.null(allNam[["PSU"]]) && engine=="survey" )  {
                         repW <- datL[!duplicated(datL[,allNam[["ID"]]]),]
                         repA <- generate.replicates(dat = repW, ID = allNam[["ID"]], wgt = allNam[["wgt"]], PSU = allNam[["PSU"]], repInd = allNam[["repInd"]], type=type )
                         rm(repW)
                     }  else  { repA <- NULL}
                  }
                  if(!is.null(allNam[["PSU"]]) || !is.null(allNam[["repWgt"]])) {doJK <- TRUE }  else {doJK <- FALSE}
    ### innere Schleife (= fuer jede Trendgruppe separat): splitten nach super splitter
                  allRes<- lapply( names(toAppl), FUN = function ( gr ) {
                      if(toCall %in% c("mean", "table"))  { allNam[["group.differences.by"]] <- attr(toAppl[[gr]], "group.differences.by") }
                      if( nchar(gr) == 0 ){ datL[,"dummyGroup"] <- "wholeGroup" ; allNam[["group"]] <- "dummyGroup" } else {allNam[["group"]] <- toAppl[[gr]] }
    ### check: Missings duerfen nur in abhaengiger Variable auftreten!          ### obere Zeile: problematisch!! "allNam" wird hier in jedem Schleifendurchlauf ueberschrieben -- nicht so superschoen!
                      noMis <- unlist ( c ( allNam[-na.omit(match(c("group", "dependent", "cross.differences"), names(allNam)))], toAppl[gr]) )
                      miss  <- which ( sapply(datL[,noMis], FUN = function (uu) {length(which(is.na(uu)))}) > 0 )
                      if(length(miss)>0) { cat(paste("Warning! Unexpected missings in variable(s) ",paste(names(miss), collapse=", "),".\n",sep=""))}
    ### check: gleichviele Imputationen je Nest und Gruppe? bei mehr als 2 gruppen zusaetzlich pruefen, ob alle paare besetzt sind (Kreuztabelle)
                      if(isTRUE(doCheck)) {                                     ### wenn check=TRUE, werden fuer jede Analyse des splitters
                         if ( length( toAppl[[gr]] ) > 1) {                     ### eine Reihe von checks durchgefuehrt
                              crsTab <- table(datL[,toAppl[[gr]]])
                              if ( length(which(crsTab < 10 )) > 0 ) {
                                   cat("Warning: small number of observations in some combinations of grouping variables:\n   Recommend to remove these group(s).\n")
                                   print(crsTab); flush.console()
                              }
                         }
                         impNes<- by(data = datL, INDICES = datL[, c(allNam[["nest"]], toAppl[[gr]]) ], FUN = function ( x ) { length(table(as.character(x[,allNam[["imp"]]])))}, simplify = FALSE)
                         laenge<- which(sapply(impNes, length) == 0)
                         if ( length(laenge ) > 0 ) {
                              cat(paste("Warning: ", length(laenge), " combination(s) of groups without any observations. Analysis most probably will crash.\n")); flush.console()
                         }
                         impNes<- table(impNes[setdiff (1:length(impNes), laenge)])
                         if(length(impNes) != 1 ) {cat("Warning: Number of imputations differ across nests and/or groups!"); print(impNes); flush.console()}
    ### check: gleichviele PSUs je nest?
                         if(!is.null(allNam[["PSU"]]))  {
                             psuNes<- table ( by(data = datL, INDICES = datL[,allNam[["nest"]]], FUN = function ( x ) { length(table(as.character(x[,allNam[["PSU"]]])))}, simplify = FALSE) )
                             if(length(psuNes) != 1 ) {cat("Warning: Number of PSUs differ across nests!"); print(psuNes)}
                         }
    ### check: sind fuer jede Gruppe alle Faktorstufen in allen nests und allen imputationen vorhanden? z.B. nicht in einer Imputation nur Jungen
                         impNes<- by(data = datL, INDICES = datL[, c(allNam[["nest"]], allNam[["imp"]]) ], FUN = checkNests, allNam=allNam, toAppl=toAppl, gr=gr, simplify = FALSE)
                         mess  <- unlist(lapply(impNes, FUN = function ( x ) { x[["mess"]]}))
                         impNes<- data.frame ( do.call("rbind", lapply(impNes, FUN = function ( x ) { unlist(lapply(x[["ret"]], FUN = length)) })) )
                         if ( !all ( sapply(impNes, FUN = function ( x ) { length(table(x)) } ) == 1) ) { cat("Warning: Number of units in at least one group differs across imputations!\n")}
    ### Achtung!! jetzt der check, der in der alten Version ueber 'checkData' gemacht wurde!
                         datL  <- do.call("rbind", by(data = datL, INDICES = datL[,c( allNam[["group"]], allNam[["nest"]], allNam[["imp"]])], FUN = checkData, allNam=allNam, toCall=toCall, separate.missing.indicator=separate.missing.indicator, na.rm=na.rm))
                         ok    <- table(datL[,"isClear"])
                         if(length(ok) > 1 ) { cat ( paste( ok[which(names(ok)=="FALSE")] , " of ", nrow(datL), " cases removed from analysis due to inconsistent data.\n",sep="")) }
                      }   else   {                                              ### hier geht die 'doCheck'-Klammer zu
                         mess <- NULL                                           ### wenn check == FALSE, wird ein leeres 'message'-Objekt zurueckgegeben
                      }
    ### nur fuer jk2.table(): "expected.values" aufbereiten ...
                      if(toCall=="table") {
                         misInd <- which(is.na(datL[,allNam[["dependent"]]]))
                         if(isTRUE(separate.missing.indicator)) {
                            if(length(misInd)>0) { datL[misInd,allNam[["dependent"]]] <- "<NA>"}
                         }  else {
                            if(length(misInd)>0) {
                               cat(paste("Warning: No seperate missing categorie was chosen. ", length(misInd), " missings were found anyhow for ",allNam[["dependent"]],". Missings will be deteted from the data.\n",sep=""))
                               if(length(misInd) == nrow(datL)) {stop()}
                               datL <- datL[-misInd,]
                            }
                         }
                         expected.values <- sort(unique(c(expected.values, names(table(datL[,allNam[["dependent"]]])))))
                      }
    ### nun wird der Datensatz zuerst nach Nests und je Nest nach Imputationen geteilt
                      if ( engine=="survey" || isFALSE(doJK)) {
                      anaA<- do.call("rbind", by(data = datL, INDICES = datL[,"isClear"], FUN = doSurveyAnalyses, allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter,
                             type=type, repA=repA, modus=modus, separate.missing.indicator = separate.missing.indicator, expected.values=expected.values, probs=probs,
                             nBoot=nBoot,bootMethod=bootMethod, formula=formula, forceSingularityTreatment=forceSingularityTreatment, glmTransformation=glmTransformation, pb=pb, toCall=toCall, doJK=doJK, poolMethod=poolMethod, useWec=useWec, refGrp=refGrp))
if(substr(as.character(datL[1,allNam[["ID"]]]),1,1 ) =="Y") {browser()}
                      }  else  {                                                ### hier muesste es jetzt ggf. nach BIFIEsurvey uebergeben werden
                      anaA<- do.call("rbind", by(data = datL, INDICES = datL[,"isClear"], FUN = doBifieAnalyses, allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter,
                             separate.missing.indicator = separate.missing.indicator, expected.values=expected.values, probs=probs, formula=formula, glmTransformation=glmTransformation, toCall=toCall))
                      }
                      if( "dummyGroup" %in% colnames(anaA) )  { anaA <- anaA[,-match("dummyGroup", colnames(anaA))] }
                      return(list(anaA=anaA, mess=mess))})
                  mess  <- unique(unlist(lapply(allRes, FUN = function ( x ) { x[["mess"]]})))
                  allRes<- do.call("rbind.fill", lapply(allRes, FUN = function ( x ) { x[["anaA"]]}))
    ### falls es messages gibt, werden die jetzt aufbereitet und ausgegeben
                  if ( length(mess)>0) {
                       nc <- nchar(mess)
                       fl <- paste(rep("#", times = nc+8), collapse="")         ### fl = first line
                       mes<- paste("###", mess, "###", sep=" ")
                       ges<- paste(fl,mes,fl,sep="\n")
                       cat("\n"); cat(ges); cat("\n")
                  }
                  rownames(allRes) <- NULL;  cat("\n")
    ### Output muss aufbereitet werden, sofern eine 'table'-Analyse ueber 'mean' gewrappt wurde, das macht die Funktion 'clearTab'
                  allRes <- clearTab(allRes, allNam = allNam, depVarOri = attr(datL, "depOri"), fc=fc, toCall=toCall, datL = datL)
                  allRes <- list(resT = list(noTrend = allRes), allNam = allNam, toCall = toCall, family=family)
                  return(allRes) }} }

compareTrends <- function ( jk2.out, only.vs.wholePop = FALSE, only.Mean.SD = only.Mean.SD ) {
          tv  <- which ( sapply (jk2.out, FUN = function ( x ) { length( grep("^trend$", x) > 0) }) > 0 )
    ### trendvariable finden
          if ( length( tv ) == 0 ) { stop ("Cannot found any trend variable.\n")}
          if ( length( tv ) > 1  ) { stop (paste ("Found more than one trend variable: '",paste(names(tv), collapse = "', '"),"'.\n",sep=""))}
          sel <- jk2.out[which(jk2.out[,names(tv)] == "trend"),]
    ### Trends nach allen Parametern getrennt vergleichen
          tr  <- do.call("rbind", by ( data = sel, INDICES = sel[,"parameter"], FUN = function ( prm ) {
                 if ( length (unique(prm[,"group"])) < 2 ) {
                      vgl <- NULL
                 }  else  {
                      spl <- data.frame ( combn(unique(prm[,"group"]),2), stringsAsFactors = FALSE)
                      vgl <- do.call("rbind", lapply ( spl, FUN = function ( gr ) {
                             prms<- prm[which(prm[,"group"] %in% gr),]
                             estD<- diff ( prms[which(prms[,"coefficient"] == "est"),"value"] )
                             seD <- sqrt ( sum(prms[which(prms[,"coefficient"] == "se"),"value"]^2) )
                             ret <- data.frame ( group = paste ("compareTrend=",gr[1],"|__|",gr[2],sep=""), prms[1:2,c("depVar", "modus", "parameter")], coefficient = c("est", "se"), value = c(estD, seD))
                             return(ret)}))
                 }
                 return(vgl) }))
    ### trend objekt gegebenenfalls reduzieren ... nur falls mit jk2.mean gecalled
          if ( isTRUE(only.Mean.SD) && unique(jk2.out[,"modus"]) == "jk2.mean" ) {
               tr <- tr[which(tr[,"parameter"] %in% c("mean", "sd")), ]
          }
          if ( isTRUE(only.vs.wholePop) ) {
               tr <- tr[grep("wholeGroup", tr[,"group"]), ]
          }
          return(tr)}

checkRegression <- function ( dat, allNam, useWec ) {
                   ch <- lapply( allNam[["independent"]], FUN = function ( i ) {
                         isKonst <- length(unique(dat[,i]))
                         if ( isKonst == 1) {
                              cat(paste("Warning: predictor '",i,"' is constant. Please check your data.\n",sep=""))
                         }
                         if ( class ( dat[,i] ) == "character" ) {
                              cat(paste("Warning: predictor '",i,"' has class 'character'. Please check your data.\n",sep=""))
                         }
                         if ( class ( dat[,i] ) %in% c("character", "factor") ) {
                              if ( isKonst > 15 && useWec == FALSE ) {
                                   cat(paste("Warning: predictor '",i,"' of class '",class ( dat[,i] ),"' has ",isKonst," levels. Please check whether this is intended.\n",sep=""))
                              }
                         } })   }                                               ### keine Rueckgabe

createLinkingError <- function  ( allNam = allNam, resT = resT, datL = datL, fc, toCall) {
          if ( is.null ( allNam[["linkErr"]] ) ) {
               cat("Note: No linking error was defined. Linking error will be defaulted to '0'.\n")
               allNam[["linkErr"]] <- "le"
               datL[,"le"]         <- 0
          }
          if ( fc == "jk2.table") {
    ### rohform des Rueckgabe-data.frames initialisieren
               le <- data.frame ( depVar = unique(resT[[1]][["out1"]][,"depVar"]), unique(datL[, c(unique(resT[[1]][["out1"]][,"depVar"]), allNam[["linkErr"]]),drop=FALSE]), stringsAsFactors = FALSE)
          }
    ### jetzt fuer "mean"
          if ( fc == "jk2.mean") {
               stopifnot (length(unique(datL[,allNam[["linkErr"]]])) == 1)
               le <- unique(datL[,allNam[["linkErr"]]])                         ### SD-difference: no linking error (=> 0!) (mit Karoline Sachse besprochen)
               le <- data.frame ( depVar = unique(resT[[1]][["out1"]][,"depVar"]), parameter = c("mean", "sd"), le = c(le,0), stringsAsFactors = FALSE)
          }
    ### jetzt fuer "glm"
          if ( fc == "jk2.glm") {
               stopifnot (length(unique(datL[,allNam[["linkErr"]]])) == 1)
               le <- unique(datL[,allNam[["linkErr"]]])
               le <- data.frame ( depVar = unique(as.character(resT[[1]][["out1"]][,"depVar"])), parameter = unique(as.character(resT[[1]][["out1"]][,"parameter"])), le = le, stringsAsFactors = FALSE)
          }
    ### jetzt fuer "quantile"
          if ( fc == "jk2.quantile") {
               stopifnot (length(unique(datL[,allNam[["linkErr"]]])) == 1)
               le <- data.frame ( unique(resT[[1]][["out1"]][,c("depVar", "parameter")]), le = unique(datL[,allNam[["linkErr"]]]), stringsAsFactors = FALSE)
          }
    ### allgemeine checks
          colnames(le)[2:3] <- c("parameter", "le")
          if ( nrow(le) > length(unique(le[,"parameter"]))) {
               stop("Linking errors must be unique for levels of dependent variable.\n")
          }
          return(le)}

conv.quantile      <- function ( dat.i , allNam, na.rm, group.delimiter, probs, nBoot,bootMethod, modus) {
                      ret  <- do.call("rbind", by(data = dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( sub.dat) {
                              if( all(sub.dat[,allNam[["wgt"]]] == 1) )  {      ### alle Gewichte sind 1 bzw. gleich
                                 ret   <- hdquantile(x = sub.dat[,allNam[["dependent"]]], se = TRUE, probs = probs,na.rm=na.rm )
                                 ret   <- data.frame (group = paste(sub.dat[1,allNam[["group"]],drop=FALSE], collapse=group.delimiter), depVar = allNam[["group"]], modus = modus, parameter = rep(names(ret),2), coefficient = rep(c("est","se"),each=length(ret)),value = c(ret,attr(ret,"se")),sub.dat[1,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                              } else {                                          ### wenn Gewichte gefordert, koennen SEs ueber Bootstrap bestimmt werden
                                 if(!is.null(nBoot)) {
                                     if(nBoot<5) {nBoot <- 5}
                                     if(bootMethod == "wQuantiles") {           ### Variante 1
                                         x     <- sub.dat[,allNam[["dependent"]]]
                                         ret   <- boot(data = x, statistic = function ( x, i) {wtd.quantile(x = x[i], weights = sub.dat[i,allNam[["wgt"]]], probs = probs,na.rm=na.rm )}, R=nBoot)
                                         ret   <- data.frame (group = paste(sub.dat[1,allNam[["group"]],drop=FALSE], collapse=group.delimiter), depVar = allNam[["group"]], modus = modus, parameter = rep(as.character(probs),2), coefficient = rep(c("est","se"),each=length(probs)), value = c(ret$t0, sapply(data.frame(ret$t), sd)), sub.dat[1,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                                     } else {                                   ### Variante 2
                                         ret   <- do.call("rbind", lapply(1:nBoot, FUN = function (b){
                                                  y   <- sample(x = sub.dat[,allNam[["dependent"]]], size = length(sub.dat[,allNam[["dependent"]]]), replace = TRUE, prob = sub.dat[,allNam[["wgt"]]]/sum(sub.dat[,allNam[["wgt"]]]))
                                                  ret <- hdquantile(x = y, se = FALSE, probs = probs,na.rm=na.rm )
                                                  return(ret)}))
                                         ret   <- data.frame (group = paste(sub.dat[1,allNam[["group"]],drop=FALSE], collapse=group.delimiter), depVar = allNam[["group"]], modus = modus, parameter = rep(as.character(probs),2), coefficient = rep(c("est","se"),each=length(probs)), value = c(wtd.quantile(x = sub.dat[,allNam[["dependent"]]], weights = sub.dat[,allNam[["wgt"]]], probs = probs,na.rm=na.rm ), sapply(data.frame(ret),sd)) , sub.dat[1,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                                     }
                                 } else {
                                     ret   <- wtd.quantile(x = sub.dat[,allNam[["dependent"]]], weights = sub.dat[,allNam[["wgt"]]], probs = probs,na.rm=na.rm )
                                     ret   <- data.frame (group = paste(sub.dat[1,allNam[["group"]],drop=FALSE], collapse=group.delimiter), depVar = allNam[["group"]], modus = modus, parameter = rep(as.character(probs),2), coefficient = rep(c("est","se"),each=length(probs)), value = c(ret, rep(NA, length(probs))) , sub.dat[1,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                                 }
                              }
                              return(ret)}))
                      ret[,"comparison"] <- NA
                      return(facToChar(ret))}


jackknife.quantile <- function ( dat.i , allNam, na.rm, type, repA, probs, group.delimiter, modus) {
                      typeS          <- recode(type, "'JK2'='JKn'")        ### typeS steht fuer type_Survey
                      design         <- svrepdesign(data = dat.i[,c(allNam[["group"]], allNam[["dependent"]]) ], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = 1, rscales = 1, repweights = repA[match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, mse = TRUE)
                      formel         <- as.formula(paste("~ ",allNam[["dependent"]], sep = "") )
                      quantile.imp   <- svyby(formula = formel, by = as.formula(paste("~", paste(allNam[["group"]], collapse = " + "))), design = design, FUN = svyquantile, quantiles = probs, return.replicates = TRUE, na.rm = na.rm)
                      molt           <- melt(data=quantile.imp, id.vars=allNam[["group"]], na.rm=TRUE)
                      molt[,"parameter"]   <- removeNonNumeric(as.character(molt[,"variable"]))
                      recString      <- paste("'",names(table(molt[,"parameter"])) , "' = '" , as.character(probs), "'" ,sep = "", collapse="; ")
                      molt[,"parameter"]   <- recode(molt[,"parameter"], recString)
                      molt[,"coefficient"] <- recode(removeNumeric(as.character(molt[,"variable"])), "'V'='est'")
                      return(facToChar(data.frame ( group = apply(molt[,allNam[["group"]],drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = allNam[["group"]], modus = modus, comparison = NA, molt[,c("parameter", "coefficient", "value", allNam[["group"]])], stringsAsFactors = FALSE))) }


conv.table      <- function ( dat.i , allNam, na.rm, group.delimiter, separate.missing.indicator , correct, expected.values, modus) {
                   table.cast <- do.call("rbind", by(data = dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( sub.dat) {
                                 prefix <- data.frame(sub.dat[1,allNam[["group"]], drop=FALSE], row.names = NULL, stringsAsFactors = FALSE )
                                 foo    <- make.indikator(variable = sub.dat[,allNam[["dependent"]]], name.var = "ind", force.indicators =expected.values, separate.missing.indikator = ifelse(separate.missing.indicator==TRUE, "always","no"))
                                 if(all(dat.i[,allNam[["wgt"]]] == 1)) {ret    <- data.frame ( prefix , desk(foo[,-1, drop = FALSE],na.rm=TRUE)[,c("Mittelwert", "std.err")], stringsAsFactors = FALSE )
                                 } else { ret    <- data.frame ( prefix , desk(foo[,-1, drop = FALSE], p.weights = sub.dat[,allNam[["wgt"]]],na.rm=TRUE)[,c("Mittelwert", "std.err")], stringsAsFactors = FALSE )}
                                 ret[,"parameter"] <- substring(rownames(ret),5)
                                 return(ret)}) )
                   if(!is.null(allNam[["group.differences.by"]]))   {
                      m            <- table.cast
                      m$comb.group <- apply(m, 1, FUN = function (ii) { crop(paste( ii[allNam[["group"]]], collapse = "."))})
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
                                      dif.iii<- data.frame(group = group, parameter = "chiSquareTest", coefficient = c("chi2","df","pValue"), value = c(chisq[["statistic"]],chisq[["parameter"]],chisq[["p.value"]]) , stringsAsFactors = FALSE )
                                      return(dif.iii)}))                        ### siehe http://www.vassarstats.net/dist2.html
                      difs[,"comparison"] <- "groupDiff"                        ### und http://onlinestatbook.com/2/tests_of_means/difference_means.html
                   }
                   ret        <- melt(table.cast, measure.vars = c("Mittelwert", "std.err"), na.rm=TRUE)
                   ret[,"coefficient"] <- recode(ret[,"variable"], "'Mittelwert'='est'; 'std.err'='se'")
                   ret        <- data.frame ( group = apply(ret[,allNam[["group"]],drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = allNam[["dependent"]], modus = modus, comparison = NA, ret[,c("coefficient", "parameter")], value = ret[,"value"], ret[,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                   if(!is.null(allNam[["group.differences.by"]]))   {return(facToChar(rbind.fill(ret,difs)))} else {return(facToChar(ret))}}


jackknife.table <- function ( dat.i , allNam, na.rm, group.delimiter, type, repA, separate.missing.indicator, expected.values, modus) {
                   dat.i[,allNam[["dependent"]]] <- factor(dat.i[,allNam[["dependent"]]], levels = expected.values)
                   typeS     <- recode(type, "'JK2'='JKn'")
                   design    <- svrepdesign(data = dat.i[,c(allNam[["group"]], allNam[["dependent"]])], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = 1, rscales = 1, repweights = repA[match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, mse = TRUE)
                   formel    <- as.formula(paste("~factor(",allNam[["dependent"]],", levels = expected.values)",sep=""))
                   means     <- svyby(formula = formel, by = as.formula(paste("~", paste(as.character(allNam[["group"]]), collapse = " + "))), design = design, FUN = svymean, deff = FALSE, return.replicates = TRUE)
                   cols      <- match(paste("factor(",allNam[["dependent"]],", levels = expected.values)",expected.values,sep=""), colnames(means))
                   colnames(means)[cols] <- paste("est",expected.values, sep="____________")
                   cols.se   <- grep("^se[[:digit:]]{1,5}$", colnames(means) )
                   stopifnot(length(cols) == length(cols.se))
                   colnames(means)[cols.se] <- paste("se____________", expected.values, sep="")
                   molt      <- melt(data=means, id.vars=allNam[["group"]], na.rm=TRUE)
                   splits    <- data.frame ( do.call("rbind", strsplit(as.character(molt[,"variable"]),"____________")), stringsAsFactors = FALSE)
                   colnames(splits) <- c("coefficient", "parameter")
                   ret       <- data.frame ( group = apply(molt[,allNam[["group"]],drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = allNam[["dependent"]], modus = modus, comparison = NA, splits, value = molt[,"value"], molt[,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE)
                   if(!is.null(allNam[["group.differences.by"]]))   {
                      m            <- ret
                      m$comb.group <- apply(m, 1, FUN = function (ii) { crop(paste( ii[allNam[["group"]]], collapse = "."))})
                      m$all.group  <- 1
                      res.group    <- tempR <- setdiff(allNam[["group"]], allNam[["group.differences.by"]])
                      if(length(res.group) == 0 ) {res.group <- "all.group"}
                      difs           <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                        if(length(tempR)>0) {
                                           datSel    <- merge(dat.i, iii[!duplicated(iii[,res.group]),res.group,drop=FALSE], by = res.group, all = FALSE)
                                        }  else  {
                                           datSel    <- dat.i
                                        }
                                        designSel <- svrepdesign(data = datSel[,c(allNam[["group.differences.by"]], allNam[["dependent"]])], weights = datSel[,allNam[["wgt"]]], type=typeS, scale = 1, rscales = 1, repweights = repA[match(datSel[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, mse = TRUE)
                                        formel    <- as.formula(paste("~", allNam[["group.differences.by"]] ,"+",allNam[["dependent"]],sep=""))
                                        tbl       <- svychisq(formula = formel, design = designSel, statistic = "Chisq")
                                        scumm     <- iii[!duplicated(iii[,res.group]),res.group,drop = FALSE]
                                        group     <- paste( paste( colnames(scumm), as.character(scumm[1,]), sep="="), sep="", collapse = ", ")
                                        dif.iii   <- data.frame(group = group, parameter = "chiSquareTest", coefficient = c("chi2","df","pValue"), value = c(tbl[["statistic"]],tbl[["parameter"]],tbl[["p.value"]]) , stringsAsFactors = FALSE )
                                        return(dif.iii)                         ### siehe http://www.vassarstats.net/dist2.html
                      } ))                                                      ### http://onlinestatbook.com/2/tests_of_means/difference_means.html
                      difs[,"comparison"] <- "groupDiff"
                   }
                   if(!is.null(allNam[["group.differences.by"]]))   {return(facToChar(rbind.fill(ret,difs)))} else {return(facToChar(ret))}}


conv.mean      <- function (dat.i , allNam, na.rm, group.delimiter, modus) {
                  deskr    <- data.frame ( do.call("rbind", by(data = dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function ( sub.dat) {
                              prefix <- sub.dat[1,allNam[["group"]], drop=FALSE]
                              if ( all(sub.dat[,allNam[["wgt"]]] == 1) )  { useWGT <- NULL}  else  {useWGT <- sub.dat[,allNam[["wgt"]]] }
                              ret    <- data.frame ( nValidUnweighted = length(na.omit(sub.dat[, allNam[["dependent"]] ])), prefix, desk(sub.dat[, allNam[["dependent"]] ], p.weights = useWGT, na.rm=na.rm)[,c("N", "N.valid", "Mittelwert", "std.err", "Varianz", "Streuung")], stringsAsFactors = FALSE)
                              names(ret) <- c( "nValidUnweighted", allNam[["group"]] , "Ncases", "NcasesValid", "mean", "se.mean", "var","sd")
                              return(ret)})), modus=modus, stringsAsFactors = FALSE)
                  if(!is.null(allNam[["group.differences.by"]]))   {
                     nCat <- table(as.character(dat.i[,allNam[["group.differences.by"]]]))
                     if ( length(nCat) < 2 ) {
                          cat(paste("Warning: Grouping variable '", allNam[["group.differences.by"]], "' only has one category within imputation and/or nest. Group differences cannot be computed. Skip computation.\n",sep=""))
                     }  else  {
                          m            <- deskr
                          m$comb.group <- apply(m, 1, FUN = function (ii) { crop(paste( ii[allNam[["group"]]], collapse = "."))})
                          m$all.group  <- 1
                          res.group    <- tempR <- setdiff(allNam[["group"]], allNam[["group.differences.by"]])
                          if(length(res.group) == 0 ) {res.group <- "all.group"}
                          kontraste    <- combn ( x = sort(unique(as.character(m[,allNam[["group.differences.by"]]]))), m = 2, simplify = FALSE)
                          difs         <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                          ret <- do.call("rbind", lapply(kontraste, FUN = function ( k ) {
                                                 if ( sum ( k %in% iii[,allNam[["group.differences.by"]]]) != length(k) ) {
                                                    cat(paste("Warning: cannot compute contrasts for 'group.differences.by = ",allNam[["group.differences.by"]],"'.\n",sep="")); flush.console()
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
                                                    dif.iii   <- data.frame(dummy, group = paste(group, paste(k, collapse = ".vs."),sep="____"), parameter = "mean", coefficient = c("est","se"), value = c(true.diff, sqrt( sum(vgl.iii[,"sd"]^2 / vgl.iii[,"nValidUnweighted"]) )) , stringsAsFactors = FALSE )
                                                    stopifnot(nrow(dif.iii)==2, nrow(vgl.iii) == 2)
                                                    dummy2    <- dif.iii[1,]
                                                    dummy2[,"coefficient"] <- "es"
                                                    dummy2[,"value"]       <- dif.iii[which(dif.iii[,"coefficient"] == "est"),"value"] / sqrt(0.5*sum(vgl.iii[,"sd"]^2))
                                                    return(dif.iii)             ### siehe http://www.vassarstats.net/dist2.html
                                                 } }))                          ### http://onlinestatbook.com/2/tests_of_means/difference_means.html
                                          return(ret)}))
                        difs[,"comparison"] <- "groupDiff"
                     }
                  }
                  deskrR   <- melt(data = deskr, id.vars = allNam[["group"]], measure.vars = setdiff(colnames(deskr), c("nValidUnweighted", "modus", allNam[["group"]]) ), na.rm=TRUE)
                  deskrR[,"coefficient"] <- recode(deskrR[,"variable"], "'se.mean'='se';else='est'")
                  deskrR[,"parameter"]   <- gsub("se.mean","mean",deskrR[,"variable"])
                  deskrR   <- data.frame ( group = apply(deskrR[,allNam[["group"]],drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = allNam[["dependent"]], modus=modus, comparison = NA, deskrR[,c( "parameter", "coefficient", "value", allNam[["group"]])], stringsAsFactors = FALSE)
                  if(!is.null(allNam[["group.differences.by"]]))   {
                      if ( length (nCat ) >1 ) {
                           return(facToChar(rbind.fill(deskrR,difs)))
                      }   else  {
                           return(facToChar(deskrR))
                      }
                  }  else {
                      return(facToChar(deskrR))
                  } }

jackknife.mean <- function (dat.i , allNam, na.rm, group.delimiter, type, repA, modus) {
          dat.i[,"N_weighted"] <- dat.i[,"N_weightedValid"] <- 1
          if( length(which(is.na(dat.i[,allNam[["dependent"]]]))) > 0 ) { dat.i[which(is.na(dat.i[,allNam[["dependent"]]])), "N_weightedValid" ] <- 0 }
          typeS<- recode(type, "'JK2'='JKn'")
          repl <- repA[ match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]]),]
          des  <- svrepdesign(data = dat.i[,c(allNam[["group"]], allNam[["dependent"]], "N_weighted", "N_weightedValid")], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = 1, rscales = 1, repweights = repl[,-1, drop = FALSE], combined.weights = TRUE, mse = TRUE)
          rets <- data.frame ( target = c("Ncases", "NcasesValid", "mean", "var"), FunctionToCall = c("svytotal","svytotal","svymean","svyvar"), formelToCall = c("paste(\"~ \", \"N_weighted\",sep=\"\")","paste(\"~ \", \"N_weightedValid\",sep=\"\")","paste(\"~ \",allNam[[\"dependent\"]], sep = \"\")","paste(\"~ \",allNam[[\"dependent\"]], sep = \"\")"), naAction = c("FALSE","TRUE","na.rm","na.rm"), stringsAsFactors = FALSE)
          ret  <- apply(rets, 1, FUN = function ( toCall ) {                    ### svyby wird dreimal aufgerufen ...
                  do   <- paste(" res <- svyby(formula = as.formula(",toCall[["formelToCall"]],"), by = as.formula(paste(\"~\", paste(allNam[[\"group\"]], collapse = \" + \"))), design = des, FUN = ",toCall[["FunctionToCall"]],",na.rm=",toCall[["naAction"]],", deff = FALSE, return.replicates = TRUE)",sep="")
                  suppressWarnings(eval(parse(text=do)))                        ### Warning erklaert in Word-Doc, wird unterdrueckt da irrelevant fuer Paket
                  resL <- melt( data = res, id.vars = allNam[["group"]], variable.name = "coefficient" , na.rm=TRUE)
                  stopifnot(length(table(resL[,"coefficient"])) == 2)
                  resL[,"coefficient"] <- recode(resL[,"coefficient"], "'se'='se'; else ='est'")
                  resL[,"parameter"]   <- toCall[["target"]]
                  attr(resL, "original") <- res
                  return(resL)})
          sds  <- do.call("rbind", by(data = dat.i, INDICES =  dat.i[,allNam[["group"]]], FUN = function (uu) {
                  namen   <- uu[1, allNam[["group"]], drop=FALSE]               ### Warning erklaert in Word-Doc, wird unterdrueckt da irrelevant fuer Paket
                  sub.rep <- repl[ match(uu[,allNam[["ID"]]], repl[,allNam[["ID"]]] ) ,  ]
                  des.uu  <- svrepdesign(data = uu[,c(allNam[["group"]], allNam[["dependent"]])], weights = uu[,allNam[["wgt"]]], type=typeS, scale = 1, rscales = 1, repweights = sub.rep[,-1, drop = FALSE], combined.weights = TRUE, mse = TRUE)
                  var.uu  <- suppressWarnings(svyvar(x = as.formula(paste("~",allNam[["dependent"]],sep="")), design = des.uu, deff = FALSE, return.replicates = TRUE, na.rm = na.rm))
                  ret     <- data.frame(namen, est = as.numeric(sqrt(coef(var.uu))), se =  as.numeric(sqrt(vcov(var.uu)/(4*coef(var.uu)))), stringsAsFactors = FALSE )
                  return(ret)}) )
          sds  <- data.frame ( melt(data = sds, id.vars = allNam[["group"]], variable.name = "coefficient" , na.rm=TRUE), parameter = "sd", stringsAsFactors = FALSE)
          resAl<- rbind(do.call("rbind",ret), sds)
          resAl<- data.frame ( group = apply(resAl[,allNam[["group"]],drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = allNam[["dependent"]], modus=modus, comparison = NA, resAl[,c("parameter","coefficient","value",allNam[["group"]])] , stringsAsFactors = FALSE)
          if(!is.null(allNam[["group.differences.by"]]))   {
             nCat <- table(as.character(dat.i[,allNam[["group.differences.by"]]]))
             if ( length(nCat) < 2 ) {
                  cat(paste("Warning: Grouping variable '", allNam[["group.differences.by"]], "' only has one category within imputation and/or nest. Group differences cannot be computed. Skip computation.\n",sep=""))
             }  else  {
                m1   <- attr(ret[[ which(rets[,"target"] == "mean") ]], "original")
                sd   <- sds[which(sds[,"coefficient"] == "est"),]
                colnames(sd) <- recode(colnames(sd), "'value'='standardabweichung'")
                m    <- merge(m1, sd, by = allNam[["group"]], all = TRUE)
                stopifnot(nrow(m) == nrow(m1))
                m$comb.group <- apply(m, 1, FUN = function (ii) {crop(paste( ii[allNam[["group"]]], collapse = "."))})
                repl1<- data.frame(t(attr(attr(ret[[which(rets[,"target"] == "mean")]], "original"), "replicates")), stringsAsFactors = FALSE )
                colnames(repl1) <- replCols <- paste("replNum", 1:ncol(repl1), sep="")
                repl1[,"comb.group"] <- rownames(repl1)
                m    <- merge(m, repl1, by = "comb.group" )
                m$all.group    <- 1
                res.group      <- tempR  <- setdiff(allNam[["group"]], allNam[["group.differences.by"]])
                if(length(res.group) == 0 ) {res.group <- "all.group"}
                kontraste      <- combn ( x = sort(unique(as.character(m[,allNam[["group.differences.by"]]]))), m = 2, simplify = FALSE)
                difs           <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                  ret <- do.call("rbind", lapply(kontraste, FUN = function ( k ) {
                                         if ( sum ( k %in% iii[,allNam[["group.differences.by"]]]) != length(k) ) {
                                              cat(paste("Warning: cannot compute contrasts for 'group.differences.by = ",allNam[["group.differences.by"]],"'.\n",sep="")); flush.console()
                                              return(NULL)                      ### Quelle fuer dieses Vorgehen:
                                         } else {                               ### Mail SW an ZKD, 07.11.2012, 17.54 Uhr, "in Absprache mit Dirk"
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
                difsL<- data.frame ( depVar = allNam[["dependent"]], melt(data = difs, measure.vars = c("dif", "se", "es") , variable.name = "coefficient" , na.rm=TRUE), modus=modus, parameter = "mean", stringsAsFactors = FALSE)
                difsL[,"coefficient"] <- recode(difsL[,"coefficient"], "'se'='se'; 'es'='es'; else = 'est'")
                difsL[,"comparison"]  <- "groupDiff"
                difsL[,"group"] <- apply(difsL[,c("group","vgl")],1,FUN = function (z) {paste(z,collapse="____")})
                resAl<- rbind(resAl,difsL[,-match("vgl", colnames(difsL))])
             }
          }
          return(facToChar(resAl)) }

giveRefgroup <- function ( refGrp) {
          if ( is.null(refGrp)) {return(".vs.wholeGroup")}
          ret <- paste0(".vs.", paste(apply(refGrp, MARGIN = 1, FUN = function ( zeile ) { paste(zeile, collapse="")}), collapse="_"))
          return(ret)}

jackknife.cov <- function (dat.i , allNam, na.rm, group.delimiter, type, repA, refGrp){
          repl<- merge(repA, dat.i, by = allNam[["ID"]], all.x = FALSE, all.y = TRUE)
          covs<- do.call("rbind", lapply(setdiff(colnames(repA),allNam[["ID"]]), FUN = function ( r ) {
                 gm <- by(data=repl[which(repl[,r]>0),], INDICES = repl[which(repl[,r]>0),allNam[["group"]]], FUN = function ( j ){wtd.mean(x=j[,allNam[["dependent"]]], weights=j[,r])})
                 ms <- as.vector(gm)
                 names(ms) <- names(gm)
                 return(data.frame ( totalMean = wtd.mean(x=repl[,allNam[["dependent"]]], weights = repl[,r] ),t(ms), stringsAsFactors=FALSE))}))
          covs<- cov(covs)[1,]
          rs  <- data.frame ( group = paste0(names(covs)[-1], giveRefgroup(refGrp)), depVar =allNam[["dependent"]], modus = NA, comparison = NA, parameter = "cov", coefficient="est", value=covs[-1], zusatz = names(covs)[-1] , stringsAsFactors = FALSE)
          colnames(rs) <- recode(colnames(rs), paste0("'zusatz'='",allNam[["group"]],"'"))
          return(rs)}

conv.cov <- function (dat.i, allNam, na.rm, group.delimiter, nBoot, refGrp){
          covs<- boot(data=dat.i, R = nBoot, statistic = function ( x, i) { c(wtd.mean(x[i,allNam[["dependent"]]], weights = x[i,allNam[["wgt"]]]), as.vector(unlist(by(data = x[i,], INDICES = x[i,allNam[["group"]]], FUN = function ( g ) { wtd.mean(g[,allNam[["dependent"]]], weights = g[,allNam[["wgt"]]])}))))})
          covs<- cov(covs$t)[1,]                                                ### Boot Method wird jetzt nicht mehr variiert, ergibt keinen Sinn
          rs  <- data.frame ( group = paste0(names(table(dat.i[,allNam[["group"]]])), giveRefgroup(refGrp)), depVar =allNam[["dependent"]], modus = NA, comparison = NA, parameter = "cov", coefficient="est", value=covs[-1], zusatz = names(table(dat.i[,allNam[["group"]]])), stringsAsFactors = FALSE)
          colnames(rs) <- recode(colnames(rs), paste0("'zusatz'='",allNam[["group"]],"'"))
          return(rs)}

### Hilfsfunktion fuer jk2.glm()
jackknife.glm <- function (dat.i , allNam, formula, forceSingularityTreatment, glmTransformation, na.rm, group.delimiter, type, repA, modus, useWec) {
if(substr(as.character(dat.i[1,allNam[["ID"]]]),1,1 ) =="Y") {browser()}
                 if ( class(glm.family) == "function" ) {
                      link <- strsplit(capture.output(str(glm.family)), split="\"")[[1]][2]
                      fam  <- recode(link, "'identity'='gaussian'; 'logit'='binomial'; 'probit'='binomial'")
                 }  else  {
                      link <- glm.family$link
                      fam  <- glm.family$family
                 }
                 sub.ana <- by(data = dat.i, INDICES = dat.i[,allNam[["group"]]], FUN = function (sub.dat) {
                            nam    <- sub.dat[1,allNam[["group"]],drop=FALSE]
                            glm.ii <- test <- glm(formula = formula, data = sub.dat, family = glm.family)
                            singular       <- names(glm.ii$coefficients)[which(is.na(glm.ii$coefficients))]
                            if(!is.null(repA)) {
                                typeS      <- recode(type, "'JK2'='JKn'")
                                design     <- svrepdesign(data = sub.dat[,c(allNam[["group"]], allNam[["independent"]], allNam[["dependent"]]) ], weights = sub.dat[,allNam[["wgt"]]], type=typeS, scale = 1, rscales = 1, repweights = repA[match(sub.dat[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, mse = TRUE)
                                if(length(singular) == 0 & forceSingularityTreatment == FALSE ) {
                                   glm.ii  <- svyglm(formula = formula, design = design, return.replicates = FALSE, family = glm.family)
                                }
                            }
                            r.squared      <- data.frame ( r.squared = var(glm.ii$fitted.values)/var(glm.ii$y) , N = nrow(sub.dat) , N.valid = length(glm.ii$fitted.values) )
                            r.nagelkerke   <- NagelkerkeR2(glm.ii)
                            summaryGlm     <- summary(glm.ii)
    ### AIC, deviance etc. nur auslesen, wenn linkfunktion NICHT 'identity'
                            if( class(glm.family) == "family" ) {
                                if (  length(grep("binomial", crop(capture.output(glm.family)))) > 0 ) {
                                      res.bl <- data.frame ( group=paste(sub.dat[1,allNam[["group"]]], collapse=group.delimiter), depVar =allNam[["dependent"]],modus = modus, parameter = c(rep(c("Ncases","Nvalid",names(glm.ii$coefficients)),2),"R2","R2nagel", "deviance", "null.deviance", "AIC", "df.residual", "df.null"),
                                                coefficient = c(rep(c("est","se"),each=2+length(names(glm.ii$coefficients))),rep("est", 7)),
                                                value=c(r.squared[["N"]],r.squared[["N.valid"]],glm.ii$coefficient,NA,NA,summaryGlm$coef[,2],r.squared[["r.squared"]],r.nagelkerke[["R2"]], test$deviance, test$null.deviance, test$aic, test$df.residual, test$df.null),sub.dat[1,allNam[["group"]], drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
                                }   else  {
                                      res.bl <- data.frame ( group=paste(sub.dat[1,allNam[["group"]]], collapse=group.delimiter), depVar =allNam[["dependent"]],modus = modus, parameter = c(rep(c("Ncases","Nvalid",names(glm.ii$coefficients)),2),"R2","R2nagel"),
                                                coefficient = c(rep(c("est","se"),each=2+length(names(glm.ii$coefficients))),rep("est", 2)),
                                                value=c(r.squared[["N"]],r.squared[["N.valid"]],glm.ii$coefficient,NA,NA,summaryGlm$coef[,2],r.squared[["r.squared"]],r.nagelkerke[["R2"]]),sub.dat[1,allNam[["group"]], drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
                                }
                            }  else  {
                                res.bl <- data.frame ( group=paste(sub.dat[1,allNam[["group"]]], collapse=group.delimiter), depVar =allNam[["dependent"]],modus = modus, parameter = c(rep(c("Ncases","Nvalid",names(glm.ii$coefficients)),2),"R2","R2nagel"),
                                          coefficient = c(rep(c("est","se"),each=2+length(names(glm.ii$coefficients))),rep("est", 2)),
                                          value=c(r.squared[["N"]],r.squared[["N.valid"]],glm.ii$coefficient,NA,NA,summaryGlm$coef[,2],r.squared[["r.squared"]],r.nagelkerke[["R2"]]),sub.dat[1,allNam[["group"]], drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
                            }
                            if(!is.null(repA) || isTRUE(useWec) ) {             ### jetzt kommt die Behandlung, wenn zusaetzlich zu JK noch singularitaet auftritt. das ueberschreibt nun die bisherige "res.bl"
                                if(length(which(is.na(glm.ii$coefficients))) > 0 ) {
                                   cat(paste("Singularity problem in regression estimation for ", length(singular)," coefficient(s): ",paste(singular, collapse = ", "),". Try workaround ... \n", sep = "")); flush.console()
                                }
                                if(isTRUE(forceSingularityTreatment) && isFALSE(useWec)) {
                                   cat("Compute coefficients in the expectation of singularities ... \n"); flush.console()
                                }
                                if(length(singular) > 0 || forceSingularityTreatment == TRUE ) {
                                   stopifnot(length(as.character(formula)) == 3 )
                                   formelNew  <- paste ( as.character(formula)[2] ," ~ ",as.character(formula)[3],sep="")
                                   if ( isFALSE(useWec) ) {
                                       cat("Unidentified bug with Nagelkerkes r^2 in singularity treatment. No r^2 is computed.\n")
                                       if ( glmTransformation == "none" )  {string <- paste("resRoh <- data.frame( withReplicates(design, quote(eatRep:::getOutputIfSingular(glm(formula = ",formelNew,", weights=.weights, family = ",fam,"(link=\"", link,"\"))))), stringsAsFactors = FALSE)",sep="")}
                                       if ( glmTransformation == "sdY" )   {string <- paste("resRoh <- data.frame( withReplicates(design, quote(eatRep:::getOutputIfSingularT1(glm(formula = ",formelNew,", weights=.weights, family = ",fam,"(link=\"", link,"\"))))), stringsAsFactors = FALSE)",sep="")}
                                       eval ( parse ( text = string ) )
                                   }  else  {                                   ### untere Zeile: Kontraste fuer jedes replicate separat bestimmen
#                                      d      <- merge(dat.i, repA, by = allNam[["ID"]], all=TRUE)
#                                      rr     <- lapply( setdiff(colnames(repA), allNam[["ID"]]), FUN = function ( reps ) {
#                                                dw <- d[which(d[,reps] > 0),]
#                                                contrasts(dw[,as.character(formula)[3]]) <- contr.wec.weighted(dw[,as.character(formula)[3]], omitted=names(table(dw[,as.character(formula)[3]]))[1], weights = dw[,reps])
#                                                eval(parse(text = paste0("re <- lm_robust(formula = as.formula(formelNew), data = dw, se_type = \"stata\", weights = ",reps,")") ) )
#                                                re <- c(re[["coefficients"]], R2 = re[["r.squared"]], Nvalid = length(re$fitted.values))
#                                                return(re)})
#                                      allReps<- do.call("rbind", rr)
#                                      thetas <- colMeans(allReps)
#                                      diffs  <- (matrix(thetas, nrow=1)[rep(1, nrow(allReps)),] - allReps)^2
#                                      resRoh <- data.frame ( theta = thetas, SE = sqrt(colSums(diffs)))
                                       contrasts(dat.i[,as.character(formula)[3]]) <- contr.wec.weighted(dat.i[,as.character(formula)[3]], omitted=names(table(dat.i[,as.character(formula)[3]]))[1], weights = dat.i[,allNam[["wgt"]]])
                                       if(!is.null(repA) ) {
                                           design1<- svrepdesign(data = dat.i[,c(allNam[["group"]], allNam[["independent"]], allNam[["dependent"]]) ], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = 1, rscales = 1, repweights = repA[match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, mse = TRUE)
                                           string1<- paste("resRoh1 <- data.frame( withReplicates(design1, quote(eatRep:::getOutputIfSingularWec(lm_robust(formula = ",formelNew,", weights=.weights, se_type =\"stata\")))), stringsAsFactors = FALSE)",sep="")
                                           eval ( parse ( text = string1 ) )
                                           contrasts(dat.i[,as.character(formula)[3]]) <- contr.wec.weighted(dat.i[,as.character(formula)[3]], omitted=names(table(dat.i[,as.character(formula)[3]]))[length(names(table(dat.i[,as.character(formula)[3]])))], weights =  dat.i[,allNam[["wgt"]]])
                                           design2<- svrepdesign(data = dat.i[,c(allNam[["group"]], allNam[["independent"]], allNam[["dependent"]]) ], weights = dat.i[,allNam[["wgt"]]], type=typeS, scale = 1, rscales = 1, repweights = repA[match(dat.i[,allNam[["ID"]]], repA[,allNam[["ID"]]] ),-1,drop = FALSE], combined.weights = TRUE, mse = TRUE)
                                           string2<- paste("resRoh2 <- data.frame( withReplicates(design2, quote(eatRep:::getOutputIfSingularWec(lm_robust(formula = ",formelNew,", weights=.weights, se_type = \"stata\")))), stringsAsFactors = FALSE)",sep="")
                                           eval ( parse ( text = string2 ) )
                                           resRoh <- rbind(resRoh1, resRoh2)[which(!duplicated(c(row.names(resRoh1), row.names(resRoh2)))),]
                                       }  else  {
                                           if (is.null(allNam[["wgt"]])) {
                                               res1   <- lm_robust(as.formula(formelNew), data = dat.i, se_type = "stata")
                                           }  else  {
                                               str1   <- paste0("res1 <- lm_robust(as.formula(formelNew), data = dat.i, weights = ",allNam[["wgt"]],", se_type = \"stata\")")
                                               eval(parse(text=str1))
                                           }
                                           contrasts(dat.i[,as.character(formula)[3]]) <- contr.wec.weighted(dat.i[,as.character(formula)[3]], omitted=names(table(dat.i[,as.character(formula)[3]]))[length(names(table(dat.i[,as.character(formula)[3]])))], weights =  dat.i[,allNam[["wgt"]]])
                                           if (is.null(allNam[["wgt"]])) {
                                               res2   <- lm_robust(as.formula(formelNew), data = dat.i, se_type = "stata")
                                           }  else  {
                                               str2   <- paste0("res2 <- lm_robust(as.formula(formelNew), data = dat.i, weights = ",allNam[["wgt"]],", se_type = \"stata\")")
                                               eval(parse(text=str2))
                                           }
                                           resRoh <- data.frame (theta = c(res1[["coefficients"]],res2[["coefficients"]], res1[["r.squared"]], length(res1[["fitted.values"]])), SE = c(res1[["std.error"]],res2[["std.error"]], NA, NA), stringsAsFactors = FALSE)
                                           rowNam <- c(names(res1[["coefficients"]]),names(res2[["coefficients"]]),"R2", "Nvalid")
                                           resRoh <- resRoh[which(!duplicated(rowNam)),]
                                           rownames(resRoh) <- rowNam[which(!duplicated(rowNam))]
                                       }
                                   }                                            ### rownames(resRoh) <- gsub("N", "Ncases", rownames(resRoh))
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

superSplitter <- function ( group=NULL, group.splits = length(group), group.differences.by = NULL, group.delimiter = "_" , dependent )  {
             group        <- as.list(group)
             names(group) <- unlist(group)
             if(max(group.splits)> length(group)) {group.splits[which(group.splits>length(group))] <- length(group)}
             group.splits <- unique(group.splits)
             superSplitti <- unlist(lapply(group.splits, FUN = function ( x ) {
                             spl <- combn(names(group),x)
                             if(class(spl) == "matrix") { spl <- as.list(data.frame(spl))} else {spl <- list(spl)}
                             spl <- unlist(lapply(spl, FUN = function ( y ) { paste(as.character(unlist(y)), collapse="________")}))
                             return(spl)}))
             superSplitti <- strsplit(superSplitti, "________")
             namen        <- unlist(lapply(superSplitti, FUN = function ( y ) { paste(y, collapse=group.delimiter)}))
             superSplitti <- lapply(superSplitti, FUN = function ( y ) {
                             if(!is.null(group.differences.by)) {if( group.differences.by %in% y ) { attr(y,"group.differences.by") <- group.differences.by}}
                             return(y)})
             names(superSplitti) <- namen
             return(superSplitti)}


dG <- function ( jk2.out , analyses = NULL, digits = 3, printDeviance ) {
    ### ggf. nach Trendgruppen getrennt
            trend <- lapply(names(jk2.out[["resT"]]), FUN = function (tr) {
                     cat(paste("       Trend group: '", tr, "'.\n",sep=""))
                     splitData <- by ( data = jk2.out[["resT"]][[tr]], INDICES = jk2.out[["resT"]][[tr]][,c("group", "depVar")], FUN = function ( spl ) {return(spl)})
			               if(is.null(analyses)) {analyses <- 1:length(splitData)}
                     for ( i in analyses) {
                           spl    <- splitData[[i]]
                           weg1   <- which ( spl[,"parameter"] %in% c("Ncases","Nvalid","R2","R2nagel", "deviance", "df.null", "df.residual", "null.deviance", "AIC"))
                           weg2   <- grep("wholePopDiff",spl[,"parameter"])
                           weg3   <- grep("^p$", spl[,"coefficient"])
                           ret    <- dcast(spl[-unique(c(weg1, weg2, weg3)),], parameter~coefficient)
                           ret[,"t.value"] <- ret[,"est"] / ret[,"se"]
                           df     <- spl[ spl[,"parameter"] == "Nvalid" & spl[,"coefficient"] == "est"  ,"value"] - nrow(ret)
                           ret[,"p.value"] <- 2*(1-pt( q = abs(ret[,"t.value"]), df = df ))
                           retNR  <- ret
                           ret    <- data.frame ( lapply(ret, FUN = function ( y ) {if(class(y)=="numeric") {y <- round(y, digits = digits)}; return(y)}), stringsAsFactors = FALSE)
                           groupNamen <- setdiff(colnames(spl), c("group","depVar","modus", "parameter", "coefficient","value", "comparison"))
                           cat ( paste( "            groups: ", paste( groupNamen, unlist(lapply(spl[1,groupNamen], as.character)), sep=" = ", collapse = "; "),"\n",sep=""))
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


getOutputIfSingularWec <- function ( glmRes ) {
                       coefs <- c(glmRes[["coefficients"]], R2 = summary(glmRes)[["r.squared"]], Nvalid = length(glmRes$fitted.values))
                       return(coefs)}

### Hilfsfunktion fuer jk2.glm() wenn Regression singulaere Terme enthaelt
### Achtung: negelkerke wird irgendwie falsch berechnet, keine Ahnung woran es liegt, wird ausgeschlossen
getOutputIfSingular <- function ( glmRes ) {
                       coefs <- na.omit(coef(glmRes))
                       rnagel<- unlist(NagelkerkeR2(glmRes))
                       names(rnagel) <- c("Nvalid", "R2nagel")
                       rnagel<- rnagel[1]
                       coefs <- c(coefs, R2 = var(glmRes$fitted.values)/var(glmRes$y), rnagel)
                       return(coefs)}

### wie oben, nur werden hier lineare Transformationen der Regressionskoeffizienten erlaubt
### unstandardisierte Koeffizienten werden genutzt, um den logit vorherzusagen
### jede person hat auf dieser logitvariablen einen wert
### standardabweichung der vorhergesagten logits + (pi^2) / 3 wird genutzt, um varianz der latenten variable zu bestimmen
### unstandardisierten logits werden durch varianz der lat. variablen geteilt, um standardisierte zu erhalten
### Achtung: negelkerke wird irgendwie falsch berechnet, keine Ahnung woran es liegt, wird ausgeschlossen
getOutputIfSingularT1<- function ( glmRes) {
                       coefs <- na.omit(coef(glmRes))
                       pred  <- sd ( glmRes$linear.predictors ) +  (pi^2)/3
                       coefs <- coefs/pred
                       rnagel<- unlist(NagelkerkeR2(glmRes))
                       names(rnagel) <- c("Nvalid", "R2nagel")
                       rnagel<- rnagel[1]
                       coefs <- c(coefs, R2 = var(glmRes$fitted.values)/var(glmRes$y), rnagel)
                       return(coefs)}

clearTab <- function ( jk2.table.output, allNam , depVarOri, fc, toCall, datL) {
            if ( fc == "jk2.table" && toCall == "mean" ) {
                 stopifnot ( all(jk2.table.output[,"parameter"] %in% c("meanGroupDiff", "wholePopDiff") == FALSE))
                 jk2 <- jk2.table.output[which(jk2.table.output[,"parameter"] == "mean"),]
                 stopifnot(length(unique(jk2[,"depVar"]))==1)
                 if(!is.null(depVarOri)) {
                     prm <- datL[which(datL[,as.character(jk2[1,"depVar"])]==1),depVarOri]
                     stopifnot(length(unique(prm))==1)
                     jk2[,"depVar"]    <- depVarOri
                 }  else  {
                     prm <- "1"
                 }
                 jk2[,"parameter"] <- unique(prm)
                 return(jk2)
            }  else  {
                 return(jk2.table.output)
            } }


### definiert die Funktion (mean, table, glm ... ) und die Methode (JK2, BRR, oder konventionell) fuer Funktion 'eatRep'
chooseFunction <- function (datI, allNam, na.rm, group.delimiter,type, repA, modus, separate.missing.indicator ,
                  expected.values, probs,nBoot,bootMethod, formula, forceSingularityTreatment, glmTransformation, pb, toCall, doJK, useWec, refGrp) {
        pb$tick(); flush.console()
        if( toCall == "mean" ) {                                                ### hier wird an die meanfunktion uebergeben
            if ( isTRUE(doJK) ) {
                 ana.i <- jackknife.mean (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, modus=modus)
            }  else  {
                 ana.i <- conv.mean (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, modus=modus)
            }
        }
        if( toCall == "table" ) {
            if ( isTRUE(doJK) ) {
                 ana.i <- jackknife.table ( dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, separate.missing.indicator = separate.missing.indicator, expected.values=expected.values, modus=modus)
            }  else  {
                 ana.i <- conv.table (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, separate.missing.indicator = separate.missing.indicator, correct=correct, expected.values=expected.values, modus=modus)
            }
        }
        if( toCall == "quantile" ) {
            if ( isTRUE(doJK) ) {
                 ana.i <- jackknife.quantile (dat.i = datI, allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, probs=probs, modus=modus)
            }  else  {
                 ana.i <- conv.quantile (dat.i = datI, allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, probs=probs, nBoot=nBoot,bootMethod=bootMethod, modus=modus)
            }
        }
        if( toCall == "cov" ) {
            if ( isTRUE(doJK) ) {
                 ana.i <- jackknife.cov (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, refGrp=refGrp)
            }  else  {
                 ana.i <- conv.cov (dat.i = datI , allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter, nBoot=nBoot, refGrp=refGrp)
            }
        }
    ### glms initialisieren, wird ueberschrieben, wenn toCall = 'glm'
        glms <- grps <- nams <- NULL
        if( toCall == "glm" ) {
    ### initiate checks specifically for regression models
            doChek<- checkRegression ( dat = datI, allNam=allNam, useWec=useWec)
            ana.i <- jackknife.glm ( dat.i = datI , allNam=allNam, formula=formula, forceSingularityTreatment=forceSingularityTreatment, glmTransformation=glmTransformation, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, modus=modus, useWec=useWec)
            glms  <- ana.i[["ori.glm"]]
            grps  <- ana.i[["nams"]]
            nams  <- lapply(grps, FUN = function ( x ) { paste(as.character(unlist(x)), collapse=group.delimiter)})
            ana.i <- ana.i[["sub.ana"]]
        }
        ana.i <- data.frame ( ana.i, datI[1,c(allNam[["nest"]], allNam[["imp"]]),drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
        return(list(ana.i=ana.i, glms=glms, nams=nams, grps=grps))}


reconstructResultsStructureGlm <- function ( group, neu, grps, group.delimiter, pooled, allNam, modus) {
        r2  <- lapply(neu[[group]], FUN = function ( x ) {var(x$fitted.values)/var(x$y)})
        r2n <- lapply(neu[[group]], FUN = function ( x ) {NagelkerkeR2(x)[["R2"]]})
        Nval<- lapply(neu[[group]], FUN = function ( x ) {length(x$fitted.values)})
        mat <- which ( unlist(lapply(grps[[1]], FUN = function ( x ) { paste(as.character(unlist(x)), collapse = group.delimiter)})) == group)
    ### achtung, bescheuert: manchmal heissen die Spalten im glm-Output 'est', manchmal 'estimate' ... bzw. manchmal 'se', manchmal 'std.error'
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
        }                                                                       ### ende Workaround
        ret <- data.frame ( group=group, depVar =allNam[["dependent"]],modus = modus, comparison = NA, parameter = c(rep(c("Ncases","Nvalid",rownames(out)),2),"R2","R2nagel"),
               coefficient = c(rep(c("est","se"),each=2+length(rownames(out))),rep("est", 2)) , value= c(NA, min(unlist(Nval)),out[,est], NA, NA, out[,se], pool.R2(unlist(r2), unlist(Nval), quiet = TRUE )[["m.pooled"]],
               pool.R2(unlist(r2n), unlist(Nval), quiet = TRUE )[["m.pooled"]]), grps[[1]][[mat]], stringsAsFactors = FALSE, row.names = NULL)
        return(ret)}

### MH wollte in grauer Vorzeit, dass unten stehendes keinen Fehler ergibt; da das aber spaeter schiefgeht, muss es hier eine Fehlermeldung geben
### 06.12.2019, SW: habe jetzt wieder Fehlermeldungen draus gemacht
checkData <- function ( sub.dat, allNam, toCall, separate.missing.indicator, na.rm) {
        if(!is.null(allNam[["PSU"]])) {
            nJkZones <- length(table(as.character(sub.dat[,allNam[["PSU"]]])))
            if(nJkZones<2)  {
               stop("Found group with less than 2 PSUs. Please check your data!\n")
               #sub.dat[,"isClear"] <- FALSE
            }
        }                                                                       ### untere Zeile: prueft; es darf GAR KEINE Missings geben
        if( (toCall == "table" & separate.missing.indicator == FALSE) | (toCall %in% c("mean", "quantile", "glm") & na.rm==FALSE ) )  {
            nObserved <- length(which(is.na(sub.dat[, allNam[["dependent"]]])))
            if(nObserved>0) {                                                   
               if ( toCall %in% c("mean", "quantile", "glm") ) {
                    stop("Found unexpected missing data in dependent variable for at least one group. Please check your data or set 'na.rm==TRUE'.\n")
                    #sub.dat[,"isClear"] <- FALSE
               }  else  {
                    cat("Warning! Found unexpected missing data in dependent variable for at least one group although 'separate.missing.indicator' was set to 'FALSE'. \n    Sure that this is intended? Try to continue execution ... \n"); flush.console()
               }
            }
        }
        if ( toCall %in% c("mean", "quantile", "glm") & na.rm==TRUE) {          ### es darf NICHT ALLES missing sein
            nMissing <- length(which(is.na(sub.dat[, allNam[["dependent"]]])))
            if(nMissing == nrow(sub.dat))  {
               stop("Some groups without any observed data. Please check your data!\n")
               #sub.dat[,"isClear"] <- FALSE
            }
        }
        return(sub.dat)}

checkNests <- function (x, allNam, toAppl, gr) {
        if(length(x[,allNam[["ID"]]]) != length(unique(x[,allNam[["ID"]]])))  {
           mess <- paste(" W A R N I N G !  '",allNam[["ID"]],"' variable is not unique within nests and imputations. Analysis will be most likely biased!",sep="")
        }   else   {
           mess <- NULL
        }
        if( length(toAppl[[gr]])>0) { ret <- lapply( toAppl[[gr]], FUN = function ( y ) {table(x[,y])}) } else {ret <- 1}
        return(list(ret=ret, mess = mess)) }


doSurveyAnalyses <- function (datL1, allNam, doJK, na.rm=na.rm, group.delimiter=group.delimiter, type=type, repA=repA, modus=modus, separate.missing.indicator = separate.missing.indicator, expected.values=expected.values,
                    probs=probs, nBoot=nBoot,bootMethod=bootMethod, formula=formula, forceSingularityTreatment=forceSingularityTreatment, glmTransformation=glmTransformation, pb=pb, toCall=toCall, poolMethod, useWec=useWec, refGrp=refGrp ) {
        if(datL1[1,"isClear"] == TRUE) {                                        ### nur fuer isClear==TRUE werden Analysen gemacht
            nrep<- table(datL1[, c(allNam[["nest"]], allNam[["imp"]])])
            nrep<- prod(dim(nrep))
    ### progress bar fuer analysen anzeigen
            cri1<- nrep > 4 & length(unique(datL1[,allNam[["ID"]]]))>2000 & doJK
            cri2<- FALSE
            if ( isFALSE(doJK) {
                 cri2<- nrep > 9 & length(unique(datL1[,allNam[["ID"]]]))>5000
            }
            if ( isTRUE(cri1) | isTRUE(cri2) ) {
                 pb  <- progress_bar$new( format = "    analyses [:bar] :percent in :elapsed", incomplete = " ", total = nrep, clear = FALSE, width= 60, show_after = 0.01)
            }  else  {
                 pb <- list()
                 pb$tick <- function (){return(NULL)}
            }
            ana <- do.call("rbind", by(data = datL1, INDICES = datL1[,allNam[["nest"]]], FUN = function ( datN ) {
                   anaI <- by(data = datN, INDICES = datN[,allNam[["imp"]]], FUN = chooseFunction, allNam=allNam, na.rm=na.rm, group.delimiter=group.delimiter,
                           type=type, repA=repA, modus=modus, separate.missing.indicator = separate.missing.indicator, expected.values=expected.values, probs=probs,
                           nBoot=nBoot,bootMethod=bootMethod, formula=formula, forceSingularityTreatment=forceSingularityTreatment, glmTransformation=glmTransformation, pb=pb, toCall=toCall, doJK, useWec=useWec, refGrp=refGrp)
                   glms <- lapply(anaI, FUN = function ( x ) { x[["glms"]]})
                   nams <- lapply(anaI, FUN = function ( x ) { x[["nams"]]})
                   grps <- lapply(anaI, FUN = function ( x ) { x[["grps"]]})
    ### hier werden jetzt, sofern nicht genestet imputiert wurde, Regressionskoeffizienten nach Methode "mice" gepoolt
    ### die Listenstruktur des Objekts muss umgedreht werden
                   if ( toCall == "glm" && poolMethod == "mice" && length(glms) > 1) {
                        aussen <- unlist(nams[[1]])
                        innen  <- names(nams)
                        for ( j in 1:length(glms)) { names(glms[[j]]) <- aussen }
                        neu    <- lapply(aussen, FUN = function ( a ) {lapply(innen, FUN = function ( i ) { glms[[i]][[a]]})})
                        pooled <- lapply(neu, FUN = function ( x ) { pool(as.mira(x)) } )
                        names(pooled) <- names(neu) <- aussen
    ### jetzt muss hier die ergebnisstruktur rekonstruiert werden und 'anaI' genannt werden
                        anaI   <- do.call("rbind", lapply(aussen, FUN = reconstructResultsStructureGlm, neu=neu, grps=grps, group.delimiter=group.delimiter, pooled=pooled, allNam=allNam, modus=modus))
                   }  else  {                                                       ### bei methode 'scalar' wird erst spaeter gepoolt
                        anaI <- do.call("rbind", lapply(anaI, FUN = function ( x ) { x[["ana.i"]]}))
                   }
                   return(anaI)}))
    ### es wird nur gepoolt, wenn es mehr als eine Imputation gibt!
            if ( poolMethod == "mice" && toCall == "glm")  {
                 retList <- ana
            }  else  {
                 if( length(table(ana[,allNam[["imp"]]])) > 1 ) {
                     retList <- jk2.pool ( datLong = ana, allNam=allNam, forceSingularityTreatment = forceSingularityTreatment, modus=modus)
                 }  else  {
                     retList <- ana[,-match(c(allNam[["nest"]], allNam[["imp"]]), colnames(ana))]
                 }
            }
    ### hier die dummy-Ergebnisstruktur erzeugen (noch nicht schoen, kann man vielleicht auch lassen ... )
        }  else  {
            retList <- NULL
        }
    ### p-Werte ergaenzen ... geschieht erst nach dem poolen der Standardfehler
        retList <- addSig (dat = retList, allNam = allNam)
        return(retList)}


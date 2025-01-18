### generiert Zeitstempel
genTS <- function() {
   timeStamp <- eatTools::removeNonNumeric(strsplit(format(Sys.time(), digits = 2L), " ")[[1]][2])
   return(timeStamp)}

### Liste mit default-Argumenten der eatRep-Hauptfunktion ... wird spaeter um das angereichert, fuer das keine defaults gesetzt sind
argl <- list(wgt = NULL, L1wgt=NULL, L2wgt=NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"), PSU = NULL, repInd = NULL, jkfac = NULL, repWgt = NULL, nest=NULL, imp=NULL,
        toCall = c("mean", "table", "quantile", "glm", "cov", "lmer", "glmer"), groups = NULL, refGrp = NULL, group.differences.by = NULL, cross.differences = FALSE, group.delimiter = "_",
        adjust=NULL, useEffectLiteR = TRUE, trend = NULL, linkErr = NULL, na.rm = FALSE, forcePooling = TRUE, boundary = 3, doCheck = TRUE, separate.missing.indicator = FALSE, expected.values = NULL, probs = NULL, nBoot = NULL, bootMethod = NULL, formula=NULL, family=NULL, formula.fixed=NULL, formula.random=NULL,
        forceSingularityTreatment = FALSE, glmTransformation = c("none", "sdY"), correct=TRUE, onlyCheck = FALSE, poolMethod = "mice", useWec = FALSE, reihenfolge = NULL, clusters=NULL, fc = NULL, isRecursive = FALSE, depOri = NULL, nCores=NULL)

### objekte auf dem NAMESPACE in Liste schreiben
captureObjectsInList <- function(env, exclude = NULL) {
            obj  <- ls(envir = env)
            if(!is.null(exclude)) {
                obj <- setdiff(obj, exclude)
            }
            a    <- list()
            for ( i in obj) {
                 x <- eval(parse(text=i), envir = env)
                 if ( is.null(x)) {
                      a[i] <- list(NULL)
                 }  else {
                      a[[i]] <- x
                 }
            }
            return(a)}
            
            
### Hilfsfunktion fuer jackknife.glm
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
### allNam <- list(wgt = "Gewichtungsvariable")
### write(createCall(hetero = TRUE, weights = "wgt", clusters = NULL, allNam=allNam), file = "c:/diskdrv/Winword/Psycho/IQB/temp/29_stan_wd/test1.r")

### Hilfsfunktion fuer eatRep: abhaengige und unabhaengige Variablen identifizieren
identify_UV_AV <- function ( a, glmerFormula)  {
          a["independent"] <- list(NULL)                                        ### initialisieren
          if(a%$$%toCall == "glm") {                                            ### fuer glm und lmer muessen abhaengge und unabhaengige Variablen aus Formel extrahiert werden
             a$dependent  <- as.character(a%$$%formula)[2]
             a$independent<- all.vars(a$formula[-2])
          }
          if(a%$$%toCall == "lmer") {
             a$independent<- unique(c(all.vars(a%$$%formula.random), all.vars(a%$$%formula.fixed)))
          }
          if(a%$$%toCall == "glmer") {
             a$dependent  <- as.character(glmerFormula)[2]
             a$independent<- all.vars(glmerFormula[-2])
          }
          return(a)}

### https://stackoverflow.com/questions/65650613/get-all-the-factor-variables-from-a-formula-call
extractFactorVarsFromFormula <- function ( formula) {
         r <- Reduce(paste0, deparse(formula))
         r2<- methods::el(regmatches(r, gregexpr("(?<=factor\\().*?(?=\\))", r, perl=TRUE)))
         return(r2)}

### Hilfsfunktion: steuert den default bei der Auswahl der sandwich-estimator
chooseSeType <- function ( se_type, clusters) {
       if ( length(se_type) > 1) {                                              ### untere zeile: defaults setzen
            if ( is.null(clusters)) { se_type <- "HC3" } else { se_type <- "CR2"}
       }
       se_type <- match.arg(arg = se_type, choices = c("HC3", "HC0", "HC1", "HC2", "CR0", "CR2"))
       return(se_type)}

### generateRandomJk1Zones (datSel, "IDSTUD", 100)
generateRandomJk1Zones <- function (datL, unit, nZones, name = "randomCluster") {
       datL  <- eatTools::makeDataFrame(datL, name = "data", minRow = 4, onlyWarn=FALSE)
       stopifnot(length(unit)==1)
       allVar<- list(ID = unit)
       allNam<- eatTools::existsBackgroundVariables(dat = datL, variable=unlist(allVar), warnIfMissing = FALSE)
       if ( "randomCluster" %in% colnames(datL)) {stop("Name '",name,"' already exists in data. Please choose an alternative name.")}
       if ( nZones >= length(unique(datL[,allNam])) ) { stop("Number of zones must not exceed number of units.")}
       if ( nZones >= length(unique(datL[,allNam])) / 5 ) {warning("Number of zones (",nZones,") is large compared to the number of distinct units (",length(unique(datL[,allNam])),").", immediate. = TRUE)}
       reps  <- length(unique(datL[,allNam])) / nZones
       zones <- rep(1:nZones, times = ceiling(reps))
       zones <- data.frame ( ID = sample(unique(datL[,allNam]), size = length(unique(datL[,allNam])), replace=FALSE), zone = zones[1:length(unique(datL[,allNam]))], stringsAsFactors = FALSE)
       colnames(zones)[2] <- name
       mdat  <- merge(datL, zones, by.x = allNam, by.y = "ID", all = TRUE)
       return(mdat)}

### Wrapper: ruft "eatRep()" mit selektiven Argumenten auf
repMean <- function(datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"), PSU = NULL, repInd = NULL, jkfac = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL,
            group.splits = length(groups), group.differences.by = NULL, cross.differences = FALSE, crossDiffSE = c("wec", "rep","old"), adjust = NULL, useEffectLiteR = FALSE, nBoot = 100,
            group.delimiter = "_", trend = NULL, linkErr = NULL, dependent, na.rm = FALSE, doCheck = TRUE, engine = c("survey", "BIFIEsurvey"), scale = 1, rscales = 1, mse=TRUE, rho=NULL, hetero=TRUE, se_type = c("HC3", "HC0", "HC1", "HC2", "CR0", "CR2"),
            clusters =NULL, crossDiffSE.engine= c("lavaan", "lm"), stochasticGroupSizes = FALSE, verbose = TRUE, progress = TRUE, nCores=NULL) {
            a    <- captureObjectsInList(env = environment(), exclude = "datL")
            ret  <- repMeanList(datL = datL, a=a)
            return(ret)}

repMeanList <- function (datL, a) {
            a$depOri <- attr(datL, "depOri")
            a$fc     <- attr(datL, "fc")
            a$crossDiffSE.engine <- match.arg(a%$$%crossDiffSE.engine, choices = c("lavaan", "lm"))
            a$cdse <- match.arg(arg = a%$$%crossDiffSE, choices = c("wec", "rep","old"))
    ### wenn adjustiert wird, koennen cross.level differenzen nur mit methode 'old' bestimmt werden
            if (!is.null(a%$$%adjust)) {
                if(is.list(a%$$%cross.differences) || isTRUE(a%$$%cross.differences)) {
                     if ( a%$$%cdse != "old") {
                         warning("To date, for adjusted means, cross-level differences can only be computed with method 'old'. Set 'crossDiffSE' to 'old'.", immediate. = TRUE)
                         a$cdse <- "old"
                     }
                }
            }
            a$type    <- car::recode(match.arg(arg = toupper(a%$$%type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'")
            a$se_type <- chooseSeType(a%$$%se_type, a%$$%clusters)
            datL      <- eatTools::makeDataFrame ( datL, name = "datL", minRow = 2, onlyWarn=FALSE)
            if ( is.null ( attr(datL, "modus"))) {
                  a$modus <- identifyMode ( name = "mean", type = car::recode(match.arg(arg = toupper(a%$$%type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'"))
            }  else  {
                  a$modus <- attr(datL, "modus")                                ### untere Zeile: damit man nicht immer copy/pasten muss, wenn man Argument uebergibt,
            }                                                                   ### also repInd=repInd, doCheck=doCheck, etc., wird das hier ueber eine eatTools-Funktion gemacht
            a$toCall <- "mean"
            ret      <- eatRep(datL = datL, a=a)
            if ( isFALSE(ret[["allNam"]][["cross.differences"]])) {return(ret)} ### wenn keine cross differences, dann wird hier gestoppt und zurueckgegeben
            if ( (is.list(a%$$%cross.differences) || a%$$%cross.differences == TRUE) && a%$$%cdse == "old" ) {
                 ret[["SE_correction"]] <- list(NULL)
                 class(ret[["SE_correction"]]) <- c("old", "list")
                 return(ret)
            }                                                                   ### hier beginnen cross-level differences
            if ( is.list(a%$$%cross.differences) || isTRUE(a%$$%cross.differences) ) {
                  toAppl<- superSplitter(group = ret[["allNam"]][["group"]], group.splits = a%$$%group.splits, group.differences.by = ret[["allNam"]][["group.differences.by"]], group.delimiter = a%$$%group.delimiter , dependent=ret[["allNam"]][["dependent"]] )
                  stopifnot(length ( toAppl ) > 1)
                  vgl <- do.call("rbind", lapply(1:length(toAppl), FUN = function ( y ) {
                         gdb <- attr(toAppl[[y]], "group.differences.by")
                         if ( is.null(gdb)) {gdb <- NA}
                         res <-  data.frame ( analysis.number = y, hierarchy.level = length(toAppl[[y]]), groups.divided.by = paste(toAppl[[y]], collapse=" + "), group.differences.by = gdb)
                         return(res)}))
                  ana <- lapply ( combinat::combn(x=vgl[,"analysis.number"], m=2, simplify=FALSE), FUN = function ( aa ) {
                         t1 <- abs(diff(vgl[ match(aa, vgl[,"analysis.number"]),"hierarchy.level"])) == 1
                         t2 <- TRUE                                             ### Vergleiche sollen nur angestellt werden, wenn Differenz der Hierarchielevel = 1 und
                         if ( vgl[min(aa),"hierarchy.level"] != 0 ) {           ### alle bis auf eine Gruppierungsvariable der hoeheren Ebene in der unteren enthalten sind, und
                              nam<- lapply(sort(aa), FUN = function ( z ) {     ### wenn Vergleiche Teilmenge dessen ist, was in 'cross.differences' angegeben ist
                                    n1 <- unlist(strsplit(as.character(vgl[z,"groups.divided.by"]), split=" |\\+"))
                                    n1 <- n1[which(nchar(n1)>0)]
                                    return(n1)})
                              if ( length(setdiff(nam[[2]], nam[[1]])) != 1) { t2 <- FALSE }
                         }
                         t3 <- TRUE
                         if (  is.list(a%$$%cross.differences) ) {
                             if ( sum(unlist(lapply(a%$$%cross.differences, FUN = function (v1) { all(sort(vgl[ match(aa, vgl[,"analysis.number"]),"hierarchy.level"]) == sort(v1))}))) == 0) {t3 <- FALSE}
                         }
                         if(isTRUE(t1) && isTRUE(t2) && isTRUE(t3)) {return(aa)} else {return(NULL)} })
                  ana <- ana[which(unlist(lapply(ana, FUN = function (l) {!is.null(l)}))==TRUE)]
                  if ( sum(abs(unlist(lapply(a%$$%cross.differences, diff))) > 1) > 0 ) {
                      warning("Computation of cross level differences using '",a%$$%cdse,"' method is only possible for differences according to adjacent levels. Non-adjacent levels will be ignored.", immediate. = TRUE)
                  }  else {
                      if ( a%$$%verbose ) {
                            if (a%$$%cdse != "wec" || is.null(ret[["allNam"]][["PSU"]])) {
                                cat(paste0("Compute cross level differences using '",a%$$%cdse,"' method.\n"))
                            }  else  {
                                cat(paste0("Compute cross level differences using '",a%$$%cdse,"' method. Assume ",car::recode(a%$$%hetero, "TRUE='heteroscedastic'; FALSE='homoscedastic'")," variances.\n"))
                            }
                      }
                  }
                  spl <- lapply(ana, FUN = function ( aa ) {
                         vgl <- vgl[aa,]
                         if (vgl[1,"hierarchy.level"] != 0) {                   ### wenn referenzdatensatz bspw. alle Maedchen sein soll, wird der jetzt nach den Gruppierungsvariablen gesplittet
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
                                     stopifnot(is.null(a%$$%wgt))               ### workaround: wenn keine gewichte spezifiziert wurden, werden sie auf 1 gesetzt und allNam$wgt ist nicht NULL
                                     message("   '",a%$$%cdse,"' method: Assume equally weighted cases.")
                                     gew <- NULL                                ### deswegen sucht er im Datensatz eine gewichtungsvariable, die nur 1en enthaelt und gibt Fehlermeldung
                                } else {                                        ### Fuer diesen Fall wird fuer 'wgt' das Argument 'NULL' uebergeben
                                     gew <- ret[["allNam"]][["wgt"]]
                                }                                               ### untere Zeilen: bloeder hotfix, notwendig wegen anderem Hotfix, weil ret$allNam$nest nicht NULL ist, selbst wenn keine nests spezifiziert werden (damit 'by' funktioniert)
                                if (is.null(a%$$%nest)) {ne <- NULL} else {ne <- ret[["allNam"]][["nest"]]}
                                if (is.null(a%$$%imp))  {im <- NULL} else {im <- ret[["allNam"]][["imp"]]}
                                if ( vgl[1,"hierarchy.level"] != 0) {           ### untere Zeile: gsub() -- in eatRep() werden '_' und Punkte aus Gruppierungsvariablen geloescht, sie muessen daher hier (VOR dem eatRep-Aufruf) bereits gegebenenfalls aus den Faktor levels geloescht werden
                                     rg <- eatTools::facToChar(d[1,nam, drop=FALSE])
                                     rg <- do.call("rbind", lapply(names(rg), FUN = function (r){data.frame ( groupName = r, groupValue = gsub("\\.", "", gsub("_", "",as.character(rg[1,r]))), stringsAsFactors = FALSE) }))
                                }  else  {
                                     rg <- NULL
                                }
                                if ( isTRUE(a%$$%hetero)) {
                                     if ( a%$$%cdse == "rep"){
                                         warning("Method 'rep' is not yet adapted for heterogeneous variances. Results might not be trustworthy.", immediate. = TRUE)
                                     }
                                }  else  {
                                     if(!is.null(ret[["allNam"]][["clusters"]])) {
                                         stop("If clusters are specified, 'hetero' must be TRUE.")
                                     }
                                }
                                if ( a%$$%cdse == "wec" ) {                     ### untere Zeile: fuer wec muss Gruppierungsvariable faktor sein, in 'eatRep()' werden Gruppierungsvariablen generell
                                     if ( !inherits(d[,grp], "factor") ) {      ### bzw. allgemein gecheckt, duerfen dort auch character sein, nur hier eben nicht, deshalb hier ein strengeres Kriterium
                                         warning("Group variable '",grp,"' must be of class 'factor' for '",a%$$%cdse,"'. Change class of '",grp,"' from '",class(d[,grp]),"' to 'factor'.", immediate. = TRUE)
                                         d[,grp] <- as.factor(d[,grp])
                                     }
                                     attr(d, "wrapperForWec") <- TRUE
                                     if ( isTRUE(a%$$%stochasticGroupSizes) && ( !is.null(ret[["allNam"]][["PSU"]]) || !is.null(ret[["allNam"]][["repWgt"]]) ) ) {
                                          message("To date, stochastic group sizes cannot be used in combination with replication mehods. Switch to fixed group sizes.")
                                          stochasticGroupSizes <- FALSE
                                     }
                                     if ( isTRUE(a%$$%stochasticGroupSizes) &&  a%$$%crossDiffSE.engine == "lm" ) {
                                          message("To date, stochastic group sizes cannot be computed with 'lm' engine. Switch to crossDiffSE.engine == 'lavaan'.")
                                          a$crossDiffSE.engine  <- "lavaan"
                                     }
                                     b <- repGlm(datL=d, ID=ret[["allNam"]][["ID"]], wgt = gew, type = a%$$%type, PSU = ret[["allNam"]][["PSU"]], repInd = ret[["allNam"]][["repInd"]],
                                                  repWgt = ret[["allNam"]][["repWgt"]], nest=ne, imp=im, trend = a%$$%trend,
                                                  formula = as.formula(paste0(ret[["allNam"]][["dependent"]] , " ~ ", grp)), doCheck = a%$$%doCheck, na.rm = a%$$%na.rm, useWec = TRUE,
                                                  scale = a%$$%scale, rscales = a%$$%rscales, mse=a%$$%mse, rho=a%$$%rho, hetero=a%$$%hetero, se_type=a%$$%se_type , crossDiffSE.engine= a%$$%crossDiffSE.engine, stochasticGroupSizes=a%$$%stochasticGroupSizes,
                                                  verbose=a%$$%verbose, progress=a%$$%progress, clusters=a%$$%clusters)
                                }  else  {                                      ### hier beginnt methode cross differences ohne wec, also 'rep' bzw. 'old'. problem: bei aufruf eatrep muss ich ja alle argumente in a$, die ich hier nicht explizit nenne, auf ihre defaults zuruecksetzen. Das bedeutet, die entsprechenden Eintraege in a$ zu loeschen
                                     g <- list(ID = ret[["allNam"]][["ID"]], wgt = gew, PSU = ret[["allNam"]][["PSU"]], repInd = ret[["allNam"]][["repInd"]],  toCall = "cov",  nest = ne, imp = im, groups = grp, refGrp =  rg, dependent = ret[["allNam"]][["dependent"]], engine = "survey", reihenfolge = ret[["allNam"]][["group"]], doCheck = FALSE, group.splits = length(grp))
                                     a <- c(a[setdiff(names(a), names(g))], g)
                                     keep <- c("ID", "wgt", "PSU", "repInd", "toCall", "nest", "imp", "groups","group.splits", "refGrp", "dependent", "engine", "reihenfolge", "type", "jkfac", "nBoot","trend", "na.rm", "modus", "scale", "rscales", "mse", "rho", "verbose", "progress", "clusters", "fc", "adjRegType", "cdse", "crossDiffSE.engine")
                                     weg  <- setdiff(names(a), keep)
                                     if(length(weg)>0) {a <- a[-eatTools::whereAre(weg, names(a), verbose=FALSE)]}
                                     b <- eatRep(datL =d, a=a)
                                }
                                b[["vgl"]]      <- vgl
                                b[["focGrp"]]   <- grp                          ### Fokus- und Referenzgruppe in Ergebnisoutput eintragen,
                                b[["refGrp"]]   <- "all"                        ### zur Weiterverarbeitung in den Reportingfunktionen
                                if(!is.null(rg)) {b[["refGrp"]] <- rg}
                         return(b)})
                  return(cld)})
                  spl <- unlist(spl, recursive=FALSE)
                  class(spl) <- c( car::recode(a%$$%cdse, "'wec'='wec_se_correction'; 'rep'='rep_se_correction'"), "list")
                  ret[["SE_correction"]] <- spl
            }
            return(ret)}

### Wrapper: ruft "eatRep()" mit selektiven Argumenten auf
repTable<- function(datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"),
            PSU = NULL, repInd = NULL, jkfac = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL, group.splits = length(groups), group.differences.by = NULL, cross.differences = FALSE, crossDiffSE = c("wec", "rep","old"),
            nBoot = 100, chiSquare = FALSE, correct = TRUE, group.delimiter = "_", trend = NULL, linkErr = NULL, dependent , separate.missing.indicator = FALSE,na.rm=FALSE, expected.values = NULL, doCheck = TRUE, forceTable = FALSE,
            engine = c("survey", "BIFIEsurvey"), scale = 1, rscales = 1, mse=TRUE, rho=NULL, verbose = TRUE, progress = TRUE, nCores=NULL ) {
            a    <- c(captureObjectsInList(env = environment(), exclude = "datL"), crossDiffSE.engine = "lavaan", adjust = list(NULL), se_type ="HC3")
            ret  <- repTableList(datL = datL, a=a)
            return(ret)}

repTableList <- function (datL, a) {
            a$crossDiffSE <- "old"                                              ### untere Zeile: warum so kompliziert? 'cross.differences' kann TRUE/FALSE oder eine Liste sein!
            if(isFALSE(a%$$%cross.differences) == FALSE) {message("To date, only method 'old' is applicable for cross level differences in frequency tables.")}
            a$modus <- identifyMode ( name = "table", type = car::recode(match.arg(arg = toupper(a%$$%type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'"))
            datL  <- eatTools::makeDataFrame ( datL, minRow = 2, onlyWarn=FALSE)
            b     <- a; b[["toCall"]] <- "table"
            b[["doCheck"]]   <- FALSE
            b[["onlyCheck"]] <- TRUE
            b[["fc"]]        <- "repTable"
            chk1  <- eatRep(datL =datL, a=b)
            if ( length(unique(datL[,chk1[["dependent"]]])) == 2 && isTRUE(all(sort(unique(datL[,chk1[["dependent"]]])) == 0:1)) && isFALSE(a%$$%forceTable)) {
                 attr(datL, "modus") <- a%$$%modus                              ### obere Zeile: wrapper! hier wird repTable ueber repMean aufgerufen; fuer eine kategorielle Variable
                 attr(datL,"fc") <- "repTable"                                  ### sind die Haeufigkeiten die Mittelwerte der Indikatoren der Faktorstufen; untere Zeilen, Achtung!! hier muessen immer zwei '&'-Zeichen gesetzt werden!!
                 b   <- a
                 for ( i in c("ID", "wgt", "PSU", "repInd", "nest", "imp", "trend", "linkErr", "dependent")) {b[[i]] <- chk1[[i]]}
                 ret <- repMeanList(datL = datL, a=b)
                 return(ret)
            }  else  {
                 if ( !is.null(a%$$%group.differences.by) && isFALSE(a%$$%chiSquare)) {
                    b   <- a
                    b[["toCall"]] <- "table"
                    b[["onlyCheck"]] <- TRUE                                    ### Funbktion wird erstmal nur zum checken benutzt
                    b[["fc"]] <- "repTable"
                    chk <- eatRep(datL=datL, a=b)
    ### missing handling muss vorneweg geschehen
                    isNa<- which ( is.na ( datL[, chk[["dependent"]] ] ))
                    if ( length ( isNa ) > 0 ) {
                         warning("Warning: Found ",length(isNa)," missing values in dependent variable '",chk[["dependent"]],"'.", immediate. = TRUE)
                         if ( isTRUE(a%$$%separate.missing.indicator) ) {
                              stopifnot ( length( intersect ( "missing" , names(table(datL[, chk[["dependent"]] ])) )) == 0 )
                              if(inherits(datL[, chk[["dependent"]] ], "factor")){# Hotfix: fuer Faktorvariablen funktioniert das einfache subsetting
                                  levOld <- levels(datL[, chk[["dependent"]] ]) ### dat[which(is.na(dat[,"var"])) ,"var"] <- "missing" nicht
                                  datL[, chk[["dependent"]] ] <- as.character(datL[, chk[["dependent"]] ])
                                  datL[isNa, chk[["dependent"]] ] <- "missing"
                                  datL[, chk[["dependent"]] ] <- factor(datL[, chk[["dependent"]] ], levels=c(levOld, "missing"))
                              }  else  {
                                  datL[isNa, chk[["dependent"]] ] <- "missing"
                              }
                         }  else  {
                              if ( isFALSE(a%$$%na.rm ) ) { stop("If no separate missing indicator is used ('separate.missing.indicator == FALSE'), 'na.rm' must be TRUE if missing values occur.\n")}
                              datL <- datL[-isNa,]
                         }
                    }                                                           ### Ende des missing handlings
                    frml<- as.formula ( paste("~ ",chk[["dependent"]]," - 1",sep="") )
                    datL[, chk[["dependent"]] ] <- as.character( datL[, chk[["dependent"]] ] )
                    matr<- data.frame ( model.matrix ( frml, data = datL) )     ### untere zeilen: Wrapper liefert objekt von repMean zurueck,
                    datL<- data.frame ( datL,  matr, stringsAsFactors = FALSE)  ### der Output muss deshalb in das Format von 'table' umgeformt werden, das macht die Funktion 'clearTab' am Ende von 'eatRep'
                    ret <- lapply ( colnames(matr), FUN = function ( dpd ) {
                           attr(datL, "modus") <- a%$$%modus
                           attr(datL,"depOri") <- chk[["dependent"]]
                           attr(datL,"fc") <- "repTable"
                           b   <- a
                           b[["dependent"]] <- dpd
                           for ( i in c("ID", "wgt", "PSU", "repInd", "nest", "imp", "trend", "linkErr")) {b[[i]] <- chk1[[i]]}
                           res <- repMeanList ( datL = datL, a = b)
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
                    a[["fc"]] <- "repTable"
                    a[["toCall"]] <- "table"
                    ret <- eatRep(datL =datL, a=a)
                 }
                 return(ret)}}


### Wrapper: ruft "eatRep()" mit selektiven Argumenten auf
repQuantile<- function(datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"),
            PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL, group.splits = length(groups),
            cross.differences = FALSE, group.delimiter = "_", trend = NULL, linkErr = NULL, dependent, probs = c(0.25, 0.50, 0.75),  na.rm = FALSE,
            nBoot = NULL, bootMethod = c("wSampling","wQuantiles") , doCheck = TRUE, 
            scale = 1, rscales = 1, mse=TRUE, rho=NULL, verbose = TRUE, progress = TRUE)  {
            modus      <- identifyMode ( name = "quantile", type = car::recode(match.arg(arg = toupper(type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'"))
            bootMethod <- match.arg ( bootMethod )                              ### repList = replacement list (welche defaultargumente ueberschrieben werden sollen)
            datL       <- eatTools::makeDataFrame ( datL, minRow = 2, onlyWarn=FALSE)
            repList   <- list(ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "quantile", engine="survey", nest = nest, imp = imp, groups = groups,
                          group.splits = group.splits, cross.differences=cross.differences, trend = trend, linkErr = linkErr, dependent = dependent, group.delimiter=group.delimiter,
                          probs=probs, na.rm=na.rm, nBoot=nBoot, bootMethod=bootMethod, doCheck=doCheck, modus=modus, scale = scale, rscales = rscales, mse=mse, rho=rho, verbose=verbose, progress=progress, clusters=NULL)
            eatRep(datL =datL, a = repList)}

### Wrapper: ruft "eatRep()" mit selektiven Argumenten auf
repGlm  <- function(datL, ID, wgt = NULL, type = c("none", "JK2", "JK1", "BRR", "Fay"),
            PSU = NULL, repInd = NULL, repWgt = NULL, nest=NULL, imp=NULL, groups = NULL,
            group.splits = length(groups), group.delimiter = "_", cross.differences = FALSE, trend = NULL, linkErr = NULL, formula,
            family=gaussian, forceSingularityTreatment = FALSE, glmTransformation = c("none", "sdY"), doCheck = TRUE, na.rm = FALSE,
            poolMethod = c("mice", "scalar") , useWec = FALSE, scale = 1, rscales = 1, mse=TRUE, rho=NULL,
            hetero=TRUE, se_type = c("HC3", "HC0", "HC1", "HC2", "CR0", "CR2"), clusters = NULL, crossDiffSE.engine= c("lavaan", "lm"), stochasticGroupSizes = FALSE, verbose = TRUE,
            progress = TRUE, nCores=NULL) {
            datL   <- eatTools::makeDataFrame ( datL, minRow = 2, onlyWarn=FALSE)
            modus  <- identifyMode ( name = "glm", type = car::recode(match.arg(arg = toupper(type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'") )
            poolMethod <- match.arg(poolMethod)
            crossDiffSE.engine <- match.arg(crossDiffSE.engine)
            se_type <- chooseSeType(se_type, clusters)                          ### repList = replacement list (welche defaultargumente ueberschrieben werden sollen)
            replList   <- list(ID=ID , wgt = wgt, type=type, PSU = PSU, repInd = repInd, repWgt = repWgt, toCall = "glm", nest = nest, imp = imp, groups = groups, group.splits = group.splits,
                          cross.differences = cross.differences, trend = trend, linkErr = linkErr, formula=formula, family=family, forceSingularityTreatment=forceSingularityTreatment, glmTransformation = glmTransformation,
                          group.delimiter=group.delimiter, na.rm=na.rm, doCheck=doCheck, modus=modus, poolMethod=poolMethod, useWec=useWec, scale = scale, rscales = rscales, mse=mse, rho=rho, hetero=hetero, se_type=se_type, crossDiffSE.engine=crossDiffSE.engine,
                          stochasticGroupSizes=stochasticGroupSizes, verbose=verbose, progress=progress, clusters=clusters, engine="survey", nCores=nCores)
            eatRep(datL =datL, a = replList)}


### Wrapper: ruft "BIFIEsurvey::BIFIE.twolevelreg()" mit selektiven Argumenten auf
repLmer  <- function(datL, ID, wgt = NULL, L1wgt=NULL, L2wgt=NULL, type = c("JK2", "JK1"),
            PSU = NULL, repInd = NULL, jkfac = NULL, rho=NULL, imp=NULL, group=NULL, trend = NULL,  dependent, formula.fixed, formula.random,
            doCheck = TRUE, na.rm = FALSE, clusters, verbose = TRUE) {
            datL   <- eatTools::makeDataFrame ( datL, minRow = 2, onlyWarn=FALSE)
            modus  <- identifyMode ( name = "lmer", type = car::recode(match.arg(arg = toupper(type), choices = c("JK2", "JK1")), "'FAY'='Fay'") )
            replList   <- list(ID=ID , wgt = wgt, L1wgt=L1wgt, L2wgt=L2wgt, type=type, PSU = PSU, repInd = repInd, jkfac = jkfac, toCall = "lmer", imp = imp, groups =group, group.splits = length(group), trend = trend, dependent=dependent,
                          formula.fixed=formula.fixed, formula.random=formula.random,engine="BIFIEsurvey", na.rm=na.rm, doCheck=doCheck, modus=modus, verbose=verbose, clusters=clusters, rho=rho)
            eatRep(datL =datL, a = replList)}

### Funktion ist nicht user-level, sondern wird von repMean, repTable, repQuantile, repGlm mit entsprechenden Argumenten aufgerufen ... a ist die Argumentenliste 'argl'
eatRep <- function (datL, a) {
          a     <- c(argl[setdiff(names(argl), names(a))], a)                   ### hier wird die von der user-level-Funktion (z.B. repMean) uebergebene Argumentenliste um das erweitert, was nicht explizit angegeben, aber von eatRep als default erwartet wird
          datL  <- eatTools::makeDataFrame(datL, name = "datL", minRow = 2, onlyWarn=FALSE)
          if ( isTRUE(a%$$%useWec) ) { a$forceSingularityTreatment <- TRUE; a$poolMethod <- "scalar"}
          if(is.null(a%$$%trend)) {a["linkErr"] <- list(NULL)}                  ### Hotfix ... sonst gibt es einen fehler, wenn kein Trend bestimmt werden soll, aber dennoch 'linkErr' spezifiziert wird
    ### rauskriegen, was der user urspruenglich aufgerufen hat: repMean, repTable, repGlm .. ? das muss nur beim ersten mal der ggf. rekursiven aufrufe geschehen
          if (is.null(a%$$%fc) && isFALSE(a%$$%onlyCheck)) {                    ### es muss nicht geschehen, wenn die Funktion nur zum checken benutzt wird
               beg   <- Sys.time()
               a$fc  <- identifyFunctionCall()                                  ### zeitschaetzung fuer CRAN/github rausnehmen
          }
          a$toCall<- match.arg(a%$$%toCall, choices = argl[["toCall"]])         ### 'oberste' Funktion suchen, die eatRep gecallt hat; zweiter Teil des Aufrufs ist dazu da, dass nicht "by" drinsteht, wenn "repMean" innerhalb einer anderen "by"-Funktion aufgerufen wird
          a$type  <- car::recode(match.arg(arg = toupper(a%$$%type), choices = c("NONE", "JK2", "JK1", "BRR", "FAY")), "'FAY'='Fay'")
          if ( a%$$%type == "NONE") {a$doJK <- FALSE }  else {a$doJK <- TRUE}
          a$engine<- match.arg(arg = a%$$%engine, choices = c("survey", "BIFIEsurvey"))
          a$glmTransformation <- match.arg(a%$$%glmTransformation, choices = argl[["glmTransformation"]])
          if(isFALSE(a%$$%forceSingularityTreatment) & a%$$%glmTransformation != "none") {
             message("'forceSingularityTreatment' was set to 'FALSE'. Please note that 'glmTransformation' is only possible if 'forceSingularityTreatment' is 'TRUE'.")
          }
          a   <- identify_UV_AV(a=a, glmerFormula=NULL)                         ### hier wird die Argumentenliste ggf. das erste Mal geaendert
          if(is.null(a%$$%groups))  {                                           ### defaults setzen, wenn gruppierungsvariable NULL war
             a$groups <- "wholeGroup"
             datL[,"wholeGroup"] <- "wholePop"
             groupWasNULL <- TRUE
          }  else  {
             groupWasNULL <- FALSE
          }
    ### wenn linkErr kein Scalar, sondern data.frame, muss es hier von 'existsBackgroundVariables' ausgenommen werden
          if ( length(a%$$%linkErr) > 1 ) {
             a$leFrame <- a%$$%linkErr
             a["linkErr"] <- list(NULL)                                         ### 'https://stackoverflow.com/questions/7944809/assigning-null-to-a-list-element-in-r'
          }  else  {
             a["leFrame"] <- list(NULL)
          }
          allVar<- c(a[c("ID", "wgt", "L1wgt", "L2wgt", "PSU", "repInd", "repWgt", "nest", "imp", "trend", "linkErr", "group.differences.by", "dependent", "independent", "adjust", "clusters")], list(group = a%$$%groups))
          allNam<- lapply(allVar, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = datL, variable=ii, warnIfMissing = TRUE, stopIfMissingOnVars = c(allVar[["PSU"]], allVar[["repInd"]]))})
          a     <- c(a[-eatTools::whereAre(names(allNam), names(a), verbose=FALSE)], allNam)
          a[["allNam"]] <- names(allNam)
    ### weil es rekursiv ist, muss der data.frame TROTZDEM in allNam angehangen werden
          if ( !is.null(a%$$%leFrame)) {a$linkErr <- a%$$%leFrame}
          if(a%$$%forceSingularityTreatment == TRUE && !is.null(a%$$%PSU) ) { a$poolMethod <- "scalar"}
    ### Konsistenz der JK-Argumente pruefen (Funktion hat keine Rueckgabe; es werden nur Fehlermeldungen ggf. ausgegeben)
          foo   <- checkJK.arguments(a)
    ### hotfix: bei wec wird ja die Gruppierungsvariable zur UV in regression ... wenn aus den Levels der Gruppierungsvariablen die '_'
    ### entfernt werden, muss das ja auch fuer die UV geschehen, sofern es sich um eine Regression fuer WEC handelt
    ### letzteres muss man erstmal rauskriegen ... geschieht hier ueber attribut des Datensatzes
          auchUV<- checkWecForUV(dat=datL, allNam = a[a%$$%allNam])
    ### check: Gruppierungsvariablen duerfen nicht numerisch sein, und personen muessen in gruppierungsvariable genestet sein
    ### da diese Funktion recht lange dauert, kann man sie weglassen, wenn in trendanalyse sich die Funktion selbst aufruft -
    ### denn der check fand ja bereits beim ersten mal statt.
          if (isFALSE(a%$$%isRecursive)) {                                      ### wird nur gemacht, wenn die Funktion sich nicht wiederholt selbst aufruft
              beg   <- Sys.time()
              datL  <- checkGroupVars ( datL = datL, allNam = a[a%$$%allNam], auchUV = auchUV)
          }
    ### check fuer adjustierungsvariablen: die duerfen nur numerisch oder dichotom sein. dasselbe gilt fuer L1- und L2-Praediktoren in multilevel regressionsmodellen mit BIFIEsurvey
    ### Achtung: ab hier wird der Datensatz in die Argumentenliste mit aufgenommen!
          a     <- checkForAdjustmentAndLmer (datL=datL, a=a, groupWasNULL=groupWasNULL)
    ### check: Manche Variablennamen sind in der Ergebnisstruktur fest vergeben und duerfen daher nutzerseitig nicht als namen von Gruppenvariablen vergeben werden
          foo   <- checkNameConvention( allNam = a[a%$$%allNam])
    ### fuer weighted effect coding muss UV ein Faktor sein, und es darf auch nur eine UV geben
          if (isTRUE(a%$$%useWec) ) {
              if ( length(a%$$%independent) != 1 ) {stop("Only one independent (grouping) variable is allowed for weighted effect coding.\n")}
              if(!inherits(a[["datL"]][,a%$$%independent],c("factor", "character", "logical", "integer"))) {stop(paste0("For weighted effect coding, independent (grouping) variable '",a%$$%independent,"' must be of class 'factor', 'character', 'logical', or 'integer'.\n"))}
          }
    ### Methode "mice" geht nur fuer nicht-genestete Imputationen
          if (a%$$%poolMethod == "mice" && !is.null(a%$$%nest) && a%$$%toCall == "glm" ) {
              message("Method 'mice' is not available for nested imputation. Switch to method 'scalar'.")
              a$poolMethod <- "scalar"
          }
    ### engine BIFIEsurvey geht noch nicht fuer alles, muss fuer manche bedingungen erstmal wieder zurueckgesetzt werden
    ### bifiesurvey verbieten, wenn genestet und/oder bei quantile + glm, und bifie survey verbieten wenn irgendwas anderes als jk2 oder jk1
    ### und bifie survey verbieten, wenn frequency table mit chi.square = TRUE und group.differences = TRUE
          a$engine<- checkEngine (a)
    ### Achtung: wenn Funktion nur zum checken genutzt wird, endet sie hier
          if ( isTRUE(a%$$%onlyCheck) ) {
              ret <- a[a%$$%allNam]
          }  else  {
    ### wie in 'defineModel': Funktion ruft sich rekursiv selber auf, wenn Trend bestimmt werden soll
              if(!is.null(a%$$%trend)) {                                        ### Achtung: hier wenn Trendberechnung geschehen soll
    ### check: wenn es trend gibt, sollte die variable mindestens 2 level haben
                  nlev<- length(unique(a$datL[,a%$$%trend]))
                  if(nlev == 1) {stop(paste0("Trend variable '",a%$$%trend,"' is constant with value '",unique(a$datL[,a%$$%trend]),"'."))}
    ### check: alle Kombinationen von faktoren in beiden datensaetzen?
                  chk5<- checkFactorLevels(a)
    ### repliziere Analyse zweimal, d.h. fuer alle Trendgruppen ... erstmal single core ... Funktion ruft sich 2x oder 3x selber auf
                  resT<- by ( data = a%$$%datL, INDICES = a$datL[,a%$$%trend], FUN = function ( subdat ) {
                         if(a%$$%verbose) {cat(paste("\nTrend group: '",subdat[1,a%$$%trend ], "'\n",sep=""))}
                         b     <- a; b["trend"] <- list(NULL)                   ### aendere die Argumentenliste
                         b$isRecursive <- TRUE                                  ### warum b <- a? Argumentliste wird nur temporaer ueberschrieben!
                         foo   <- eatRep(datL=subdat, a=b)
                         foo[["resT"]][["noTrend"]][,a%$$%trend] <- subdat[1,a%$$%trend ]
                         foo[["allNam"]][["trend"]] <- a%$$%trend
                         return(list(out1=foo[["resT"]][["noTrend"]], out2=foo[["allNam"]]))}, simplify = FALSE)
    ### Trends wurden berechnet: 'resT' aufbereiten, Linkingfehler ranmatchen ... wenn es keinen linkingfehler gibt, wird er auf 0 gesetzt
    ### Hotfix (war in alter eatRep-Version auch schon so): rekursiver Aufruf bei Trends setzt fuer jedes Jahr trend = NULL und linkErr = NULL
    ### Linkerror muss daher jetzt wieder rekonstruiert werden
                  a      <- keepNonNULL(list1 = a, list2 = resT[[1]][["out2"]])
                  a$allNam <- c(names(resT[[1]][["out2"]]), "cross.differences")
                  a$le     <- createLinkingError ( resT = resT, a=a)
                  out3   <- lapply(resT, FUN = function ( k ) { k[["out1"]]})
                  ret    <- list(resT = out3, allNam = a[a%$$%allNam], toCall = a%$$%toCall, family=a%$$%family, le=a%$$%le)
                  return(ret)
              }  else {
    ### obere Zeile: Ende der inneren Schleife (= Ende des Selbstaufrufs wegen trend) ... das untere wird nun fuer jeden Aufruf abgearbeitet
                  if( length( setdiff ( a%$$%group.differences.by,a%$$%group)) != 0) {stop("Variable in 'group.differences.by' must be included in 'groups'.\n")}
    ### Anzahl der Analysen aufgrund mehrerer Hierarchieebenen definieren ueber den 'super splitter' und Analysen einzeln (ueber 'lapply') starten
                  toAppl<- superSplitter(group = a%$$%group, group.splits = a%$$%group.splits, group.differences.by = a%$$%group.differences.by, group.delimiter = a%$$%group.delimiter , dependent=a%$$%dependent )
                  if(a%$$%verbose){cat(paste(length(toAppl)," analyse(s) overall according to: 'group.splits = ",paste(a%$$%group.splits, collapse = " ") ,"'.", sep=""))}
    ### bei mehr als einer Analyse: genauere Ausgabe auf Konsole printen und consistency checks
                  ret   <- createAnalysisInfTable(toAppl=toAppl, verbose=a%$$%verbose, allNam=a[a%$$%allNam])
    ### cross differences defaulten wenn vom user nicht definiert
                  a     <- setCrossDifferences (a=a)
    ### Achtung: wenn keine Gruppen und/oder Nests und/oder Imputationen spezifiziert sind, erzeuge Variablen mit Werten gleich 1, damit by() funktioniert!
                  beg   <- Sys.time()                                           ### untere Funktion veraendert ggf. das Datensatzobjekt und allNam
                  a     <- createLoopStructure(a=a)
    ### check: wenn cross.differences gemacht werden sollen, dann muessen die faktor levels aller gruppierungsvariablen disjunkt sein (siehe Mail Benjamin, 13.11.2019, 18.11 Uhr)
                  if(!is.null(a%$$%cross.differences)) {
                      if(length(a%$$%group)>1) {
                         lev <- unlist(lapply(a%$$%group, FUN = function ( v ) { unique(as.character(a$datL[,v]))}))
                         if (length(lev) != length(unique(lev))) {stop("Factor levels of grouping variables are not disjunct.\n")}
                      }
                  }
    ### check: abhaengige Var. numerisch?
                  if(a%$$%toCall %in% c("mean", "quantile", "glm")) {
                     if(!inherits(a[["datL"]][,a%$$%dependent],  c("integer", "numeric"))) {
                         stop(paste0("Dependent variable '",a%$$%dependent,"' has to be of class 'integer' or 'numeric'.\n"))
                     }
                  }
    ### replicates zuweisen bzw. erzeugen
                  a$repA  <- assignReplicates ( a=a )
    ### innere Schleife (= fuer jede Hierachieebene separat): splitten nach super splitter.
                  allRes<- innerLoop(toAppl=toAppl, ret=ret, a=a)
                  if(a%$$%verbose){cat("\n")}
    ### Output muss aufbereitet werden, sofern eine 'table'-Analyse ueber 'mean' gewrappt wurde, das macht die Funktion 'clearTab'
                  allRes <- clearTab(allRes, allNam = a[a%$$%allNam], depVarOri = a%$$%depOri, fc=a%$$%fc, toCall=a%$$%toCall, datL = a%$$%datL)
                  allRes <- prepForReport2(out=allRes, info = ret, allNam = a[a%$$%allNam])
                  allRes <- list(resT = list(noTrend = allRes), allNam = a[a%$$%allNam], toCall = a%$$%toCall, family=a%$$%family)
                  return(allRes) }} }

### Hilfsfunktion fuer eatRep(): bereitet output fuer report2() vor. Zum einen werden die Analyse-IDs ergaenzt, und es werden, falls vorhanden, hierarchieebenen ergaenzt
prepForReport2 <- function(out, info, allNam) {
       if (!"comparison" %in% colnames(out)) {
           out[,"comparison"] <- "none"
       }  else  {
           out[,"comparison"] <- car::recode(out[,"comparison"], "NA='none'")
       }
       if(!is.null(allNam[["group"]]) && !identical (allNam[["group"]], "wholeGroup") ) { for ( g in allNam[["group"]]) {out[,g] <- car::recode(out[,g], "NA='total'")  }}
       timeStamp <- genTS()                                                     ### erstmal IDs fuer Punktschaetzer (ohne vergleiche)
    ### erstmal IDs fuer Punktschaetzer (ohne vergleiche)
       out[,"row"] <- 1:nrow(out)
       if(!"none" %in% out[,"comparison"]) {                                    ### hotfix: bei WEC oder rep methode fuer cross-level diffs kann in ergebnisstruktur ggf. kein Eintrag "none" stehen, dann soll alles genommen werden.
           stopifnot(length(unique(out[,"comparison"])) ==1 )
           point  <- out
       }  else  {
           point  <- out[which(out[,"comparison"] == "none"),]
       }
       point     <- data.frame ( type = "point", do.call("rbind", by(data = point, INDICES = point[,c("group", "depVar")], FUN = function (sg) {
                    stopifnot(nrow(unique(sg[,c("parameter", "coefficient")])) == nrow(sg))
                    sg[,"id"] <- paste(timeStamp, min(sg[,"row"]), sep="_")
                    return(sg)})), stringsAsFactors = FALSE)
    ### falls es groupdifferences gibt, jetzt dafuer ids vergeben
       gd        <- out[which(out[,"comparison"] != "none"),]
       if ( nrow(gd)>0 && "none" %in% out[,"comparison"]) {
            gd   <- data.frame ( type = "point", do.call("rbind", by(data = gd, INDICES = gd[,c("group", "depVar")], FUN = function (sg) {
                    suppressWarnings(sg1 <- data.frame ( unique(sg[, c("parameter", setdiff(allNam[["group"]], allNam[["group.differences.by"]])), drop=FALSE]), strsplit(sg[1,allNam[["group.differences.by"]]], ".vs.| - ")[[1]], stringsAsFactors = FALSE))
                    colnames(sg1)[ncol(sg1)] <- allNam[["group.differences.by"]]
                    sg2 <- merge(sg1, point, by = colnames(sg1), all=FALSE)
                    stopifnot(length(unique(sg2[,"id"])) ==2 || unique(sg1[,"parameter"]) == "chiSquareTest")
                    sg  <- data.frame ( sg, id = paste(timeStamp, min(sg[,"row"]), sep="_"), unit_1 = sort(unique(sg2[,"id"]))[1], unit_2 = sort(unique(sg2[,"id"]))[2], stringsAsFactors = FALSE)
                    return(sg) })), stringsAsFactors = FALSE)
            point<- plyr::rbind.fill(point, gd)
       }
    ### Hierarchieebenen ergaenzen ... mergeAttr kann spaeter ohne warnungen, wenn sich das als ok erwiesen hat
       if(!is.null(info)) {
            inf  <- prepareInfo(info, point, allNam)
            point<- merge(point, inf, by = allNam[["group"]], all.x=TRUE, all.y = FALSE)
       }
       return(point)}

prepareInfo <- function(info, point, allNam) {
       vals <- info[,"groups.divided.by"]
       vars <- unique(unlist(strsplit(vals, " \\+ ")))
       inf  <- do.call("rbind", by(data = info, INDICES = info[,"analysis.number"], FUN = function(z) {
               l1 <- list()
               for ( v in vars) {
                   if ( grepl(v, z[["groups.divided.by"]])) {
                       levs <- setdiff(unique(point[which(point[,"comparison"] == "none"),v]), "total")
                       if ( !is.null(allNam[["group.differences.by"]]) && v == allNam[["group.differences.by"]]) {levs <- c(levs, unlist(lapply(combinat::combn(x=levs, m=2, simplify=FALSE), FUN = function (xx) {c(paste(rev(xx), collapse=" - "), paste(xx, collapse=" - "))})))}
                       l1[[v]] <- levs
                   } else {
                       l1[[v]] <- "total"
                   }
               }
               l2 <- expand.grid(l1, stringsAsFactors = FALSE)
               suppressWarnings(l2 <- data.frame ( l2, z[,"hierarchy.level", drop=FALSE], stringsAsFactors = FALSE))
               return(l2)}))
       return(inf)}

### Hilfsfunktion fuer eatRep(): check whether factor levels are equal between trend groups. Funktion hat keine Rueckgabe, gibt nur ggf. Fehlermeldung
checkFactorLevels <- function(a) {
       for ( i in c("group", "trend", "datL")) { assign(i, a[[i]]) }
       if (!is.null(group)) {
             foo <- lapply(group, FUN = function ( gr ) {
                    ch <- by(data = datL, INDICES = datL[,trend], FUN = function ( subdat ) { table(subdat[,gr]) }, simplify = FALSE )
                    cmb<- combinat::combn(x=1:length(ch), m=2, simplify=FALSE)
                    ch1<- all(unlist(lapply(cmb,FUN = function (y) {isTRUE(all.equal(sort(names(ch[[y[1]]])), sort(names(ch[[y[2]]]))))})))
                    if(!ch1) {
                        tab <- table(datL[,c(gr,trend)])
                        warning(paste0("Levels of grouping variable '",gr,"' do not match between trend groups: \n", eatTools::print_and_capture(tab, 5) ), immediate. = TRUE)
                    }})}}

### Hilfsfunktion fuer eatRep()
createAnalysisInfTable <- function(toAppl, verbose, allNam) {
         if ( length ( toAppl ) > 1) {
               ret <- do.call("rbind", lapply(1:length(toAppl), FUN = function ( y ) {
               gdb <- attr(toAppl[[y]], "group.differences.by")
               if ( is.null(gdb)) {gdb <- NA}
               res <-  data.frame ( analysis.number = y, hierarchy.level = length(toAppl[[y]]), groups.divided.by = paste(toAppl[[y]], collapse=" + "), group.differences.by = gdb)
    ### Achtung!! adjustiert werden kann nicht fuer die oberste Hierarchieebene (hierarchy level 0) ... fuer analysen auf dieser Ebene wird adjust zu NULL
               if ( !is.null(allNam[["adjust"]])) {
                      res[,"adjust"] <- car::recode(res[,"hierarchy.level"], "0='FALSE'; else = 'TRUE'")
               }
               return(res)}))
               if(verbose){cat("\n \n"); print(ret, row.names=FALSE)}
         }  else  {ret <- NULL}
         return(ret)}

### Hilfsfunktion fuer eatRep()
innerLoop <- function (toAppl, ret, a=a)  {
       allRes<- do.call(plyr::rbind.fill, lapply( names(toAppl), FUN = function ( gr ) {
                a[["str1"]]  <- createInfoString (ai = ret, toCall=a%$$%toCall, gr=gr, toAppl=toAppl)
                if(a%$$%toCall %in% c("mean", "table"))  {                      ### obere Zeile: string snippet fuer fortschrittsanzeige kreieren
                   if ( is.null(attr(toAppl[[gr]], "group.differences.by"))) {
                        a["group.differences.by"] <- list(NULL)
                   }  else  {
                        a$group.differences.by <- attr(toAppl[[gr]], "group.differences.by")
                   }
                }
                if( nchar(gr) == 0 ){
                    a$datL[,"dummyGroup"] <- "wholeGroup"
                    a$group     <- "dummyGroup"                                 ### problematisch!! "allNam" wird hier in jedem Schleifendurchlauf ueberschrieben
                    a["adjust"] <- list(NULL)                                   ### auf der obersten Hierarchieebene findet keine Adjustierung statt
                } else {
                    a$group <- toAppl[[gr]]
                }
    ### check: Missings duerfen nur in abhaengiger Variable auftreten!
                noMis <- unlist ( c ( a[a%$$%allNam][-na.omit(match(c("group", "dependent", "cross.differences"), names(a[a%$$%allNam])))], toAppl[gr]) )
                miss  <- which ( sapply(a[["datL"]][,noMis], FUN = function (uu) {length(which(is.na(uu)))}) > 0 )
                if(length(miss)>0) { warning("Unexpected missings in variable(s) ",paste(names(miss), collapse=", "),".", immediate. = TRUE)}
    ### check: gleichviele Imputationen je Nest und Gruppe? bei mehr als 2 gruppen zusaetzlich pruefen, ob alle paare besetzt sind (Kreuztabelle)
                beg   <- Sys.time()
                chk3  <- checkImpNest(toAppl = toAppl, gr=gr, a=a)
    ### nur fuer repTable(): "expected.values" aufbereiten ... Funktion veraendert Datensatz und expected.values
                a     <- prepExpecVal (a=a)
                b     <- a[-match("datL", names(a))]                            ### Datensatz aus Argumentliste entfernen
    ### nun wird der Datensatz zuerst nach Nests und je Nest nach Imputationen geteilt
                if ( a%$$%engine=="survey" || isFALSE(a%$$%doJK)) {
                     anaA<- doSurveyAnalyses(datL1 = a[["datL"]], a = b)
                }  else  {                                                      ### hier muesste es jetzt ggf. nach BIFIEsurvey uebergeben werden
                     anaA<- doBifieAnalyses(dat.i = a[["datL"]], a = b)
                }
                if( "dummyGroup" %in% colnames(anaA) )  { anaA <- anaA[,-match("dummyGroup", colnames(anaA))] }
                return(anaA)}))
       rownames(allRes) <- NULL
       return(allRes)}

### Hilfsfunktion fuer eatRep()
identifyFunctionCall <- function() {
          i     <- 0                                                            ### wenn wec genutzt werden soll, geht das nur ueber den Hotfix, der auch fuer Singularitaet verwendet wird
          fc    <- NULL                                                         ### initialisieren
          while ( !is.null(sys.call(i))) { fc <- c(fc, eatTools::crop(unlist(strsplit(deparse(sys.call(i))[1], split = "\\("))[1])); i <- i-1  }
    ### unten stehende Zeile deshalb, damit der Nutzer sowohl repMean(...) als auch eatRep::repMean(...) als auch eatRep:::repMean(...) aufrufen kann
          fc   <- unlist(lapply(strsplit(fc, ":"), FUN = function ( l ) {l[length(l)]}))
          fc   <- fc[max(which(fc %in% c("repMean", "repTable", "repGlm", "repQuantile", "repLmer", "repGlmer")))]
          return(fc)}

checkRegression <- function ( dat, allNam, useWec ) {
                   ch <- lapply( allNam[["independent"]], FUN = function ( i ) {
                         isKonst <- length(unique(dat[,i]))
                         if ( isKonst == 1) {
                              warning("Predictor '",i,"' is constant. Please check your data.", immediate. = TRUE)
                         }
                         if ( inherits ( dat[,i],  "character" )) {
                              warning("Predictor '",i,"' has class 'character'. Please check your data.", immediate. = TRUE)
                         }
                         if ( inherits ( dat[,i],  c("character", "factor") )) {
                              if ( isKonst > 15 && isFALSE(useWec) ) {
                                   warning("Predictor '",i,"' of class '",class ( dat[,i] ),"' has ",isKonst," levels. Please check whether this is intended.", immediate. = TRUE)
                              }
                         } })   }                                               ### keine Rueckgabe

createLinkingError <- function  ( resT = resT, a=a) {
          for ( i in names(a)) { assign(i, a[[i]]) }
    ### erstmal checken, ob linkingfehler als scalar oder als data.frame uebergeben wurden ... im zweiten Fall finden die checks erst in computeTrend() statt!
          if (length(linkErr) > 1) {                                            ### data.frame: objekt wird durchgeschleift, checks finden spaeter statt!
              le <- linkErr
              attr(le, "linkingErrorFrame") <- TRUE                             ### das wird gemacht, damit die check-Funktion in computeTrend() erkennt, ob es urspruenglich ein data.frame war!
              return(le)
          }  else  {
    ### hier beginnt scalar methode: nur fuer zwei MZP geeignet!
              if ( is.null ( linkErr ) ) {
                   message("Note: No linking error was defined. Linking error will be defaulted to '0'.")
                   linkErr     <- "le"
                   datL[,"le"] <- 0
              }
    ### Achtung! Wenn Linkingfehler als Spalte im datensatz uebergeben, impliziert das, das die Linkingfehler fuer alle Paare von zeitpunkten
    ### gleich sind, also LE fuer 2011 vs. 2016 == 2011 vs. 2021 == 2016 vs. 2021 ... deshalb hier als Schleife ueber alle Paare!!
    ### rohform des Rueckgabe-data.frames initialisieren
              times<- sort(unique(datL[,trend]))
              paare<- combinat::combn(x=times, m=2, simplify=FALSE)
              les  <- do.call("rbind", lapply(paare, FUN = function (p ) {
                      init <- data.frame ( trendVar = trend, trendLevel1 = p[1] , trendLevel2 = p[2], stringsAsFactors = FALSE)
                      if ( fc == "repTable") {
                           le <- data.frame ( init, depVar = unique(resT[[1]][["out1"]][,"depVar"]), unique(datL[, c(unique(resT[[1]][["out1"]][,"depVar"]), linkErr),drop=FALSE]), stringsAsFactors = FALSE)
                      }
    ### jetzt fuer "mean"
                      if ( fc == "repMean") {
                           stopifnot (length(unique(datL[,linkErr])) == 1)
                           le <- unique(datL[,linkErr])                         ### SD-difference: no linking error (=> 0!) (mit Karoline Sachse besprochen)
                           le <- data.frame ( init, depVar = unique(resT[[1]][["out1"]][,"depVar"]), parameter = c("mean", "sd"), le = c(le,0), stringsAsFactors = FALSE)
                      }
    ### jetzt fuer "glm"
                      if ( fc %in% c("repGlm", "repLmer", "repGlmer")) {
                           stopifnot (length(unique(datL[,linkErr])) == 1)
                           le <- unique(datL[,linkErr])
                           le <- data.frame ( init, depVar = unique(as.character(resT[[1]][["out1"]][,"depVar"])), parameter = unique(as.character(resT[[1]][["out1"]][,"parameter"])), le = le, stringsAsFactors = FALSE)
                      }
    ### jetzt fuer "quantile"
                      if ( fc == "repQuantile") {
                           stopifnot (length(unique(datL[,linkErr])) == 1)
                           le <- data.frame ( init, unique(resT[[1]][["out1"]][,c("depVar", "parameter")]), le = unique(datL[,linkErr]), stringsAsFactors = FALSE)
                      }
    ### allgemeine checks
                      colnames(le)[5:6] <- c("parameter", "le")
                      if ( nrow(le) > length(unique(le[,"parameter"]))) {
                           stop("Linking errors must be unique for levels of dependent variable.\n")
                      }
                      return(le)}))
          }
          return(les)}
          

conv.quantile      <- function ( dat.i , a) {
                      for ( i in names(a)) { assign(i, a[[i]]) }
                      ret  <- do.call("rbind", by(data = dat.i, INDICES = dat.i[,group], FUN = function ( sub.dat) {
                              if( all(sub.dat[,wgt] == 1) )  {                  ### alle Gewichte sind 1 bzw. gleich
                                 ret   <- Hmisc::hdquantile(x = sub.dat[,dependent], se = TRUE, probs = probs,na.rm=na.rm )
                                 ret   <- data.frame (group = paste(sub.dat[1,group,drop=FALSE], collapse=group.delimiter), depVar = dependent, modus = modus, parameter = rep(names(ret),2), coefficient = rep(c("est","se"),each=length(ret)),value = c(ret,attr(ret,"se")),sub.dat[1,group,drop=FALSE], stringsAsFactors = FALSE)
                              } else {                                          ### wenn Gewichte gefordert, koennen SEs ueber Bootstrap bestimmt werden
                                 if(!is.null(nBoot)) {
                                     if(nBoot<5) {nBoot <- 5}
                                     if(bootMethod == "wQuantiles") {           ### Variante 1
                                         x     <- sub.dat[,dependent]
                                         ret   <- boot::boot(data = x, statistic = function ( x, i) {Hmisc::wtd.quantile(x = x[i], weights = sub.dat[i,wgt], probs = probs,na.rm=na.rm )}, R=nBoot)
                                         ret   <- data.frame (group = paste(sub.dat[1,group,drop=FALSE], collapse=group.delimiter), depVar = dependent, modus = modus, parameter = rep(as.character(probs),2), coefficient = rep(c("est","se"),each=length(probs)), value = c(ret$t0, sapply(data.frame(ret$t), sd)), sub.dat[1,group,drop=FALSE], stringsAsFactors = FALSE)
                                     } else {                                   ### Variante 2
                                         ret   <- do.call("rbind", lapply(1:nBoot, FUN = function (b){
                                                  y   <- sample(x = sub.dat[,dependent], size = length(sub.dat[,dependent]), replace = TRUE, prob = sub.dat[,wgt]/sum(sub.dat[,wgt]))
                                                  ret <- Hmisc::hdquantile(x = y, se = FALSE, probs = probs,na.rm=na.rm )
                                                  return(ret)}))
                                         ret   <- data.frame (group = paste(sub.dat[1,group,drop=FALSE], collapse=group.delimiter), depVar = dependent, modus = modus, parameter = rep(as.character(probs),2), coefficient = rep(c("est","se"),each=length(probs)), value = c(Hmisc::wtd.quantile(x = sub.dat[,dependent], weights = sub.dat[,wgt], probs = probs,na.rm=na.rm ), sapply(data.frame(ret),sd)) , sub.dat[1,group,drop=FALSE], stringsAsFactors = FALSE)
                                     }
                                 } else {
                                     ret   <- Hmisc::wtd.quantile(x = sub.dat[,dependent], weights = sub.dat[,wgt], probs = probs,na.rm=na.rm )
                                     ret   <- data.frame (group = paste(sub.dat[1,group,drop=FALSE], collapse=group.delimiter), depVar = dependent, modus = modus, parameter = rep(as.character(probs),2), coefficient = rep(c("est","se"),each=length(probs)), value = c(ret, rep(NA, length(probs))) , sub.dat[1,group,drop=FALSE], stringsAsFactors = FALSE)
                                 }
                              }
                              return(ret)}))
                      ret[,"comparison"] <- NA
                      return(eatTools::facToChar(ret))}


jackknife.quantile <- function ( dat.i , a) {
     for ( i in names(a)) { assign(i, a[[i]]) }
     typeS   <- car::recode(type, "'JK2'='JKn'")                                ### typeS steht fuer type_Survey
     design  <- survey::svrepdesign(data = dat.i[,c(group, dependent) ], weights = dat.i[,wgt], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(dat.i[,ID], repA[,ID] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
     formel  <- as.formula(paste("~ ",dependent, sep = "") )
    ### Hotfix: return.replicates = FALSE gesetzt, weil es sonst ab survey version 4.1-1 eine fehlermeldung gibt ... weiss nicht, ob man die replicates spaeter nochmal braucht, ich glaube nicht
     quant   <- survey::svyby(formula = formel, by = as.formula(paste("~", paste(group, collapse = " + "))), design = design, FUN = survey::svyquantile, quantiles = probs, return.replicates = FALSE, na.rm = na.rm)
     molt    <- eatTools::facToChar(reshape2::melt(data=quant, id.vars=group, na.rm=FALSE))
     molt[,"parameter"]  <- eatTools::crop(eatTools::removePattern(eatTools::removePattern(molt[,"variable"],paste0("se.", dependent)),dependent), char=".")
     molt    <- do.call("rbind", plyr::alply(molt, .margins = 1, .fun = function (zeile) {
                coef <- eatTools::removePattern(eatTools::crop(eatTools::removePattern(zeile[["variable"]], zeile[["parameter"]]), "."), dependent)
                if ( coef=="") {coef <- "est"} else {stopifnot(coef == "se."); coef <- "se"}
                zeile[["coefficient"]] <- coef
                return(zeile)}))
     return(eatTools::facToChar(data.frame ( group = apply(molt[,group,drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = dependent, modus = paste(modus,"survey", sep="__"), comparison = NA, molt[,c("parameter", "coefficient", "value", group)], stringsAsFactors = FALSE))) }


conv.table      <- function ( dat.i , a) {
     for ( i in names(a)) { assign(i, a[[i]]) }
     tabs <- do.call("rbind", by(data = dat.i, INDICES = dat.i[,group], FUN = function ( sub.dat) {
             prefix <- data.frame(sub.dat[1,group, drop=FALSE], row.names = NULL, stringsAsFactors = FALSE )
             foo    <- make.indikator(variable = sub.dat[,dependent], name.var = "ind", force.indicators =expected.values, separate.missing.indikator = "no")
             if (all(dat.i[,wgt] == 1)) {wgts <- NULL } else { wgts <- sub.dat[,wgt]}
             ret    <- data.frame ( prefix , eatTools::descr(foo[,-1, drop = FALSE],p.weights = wgts, na.rm=TRUE)[,c("Mean", "std.err")], stringsAsFactors = FALSE )
             ret[,"parameter"] <- substring(rownames(ret),5)
             return(ret)}) )
     Ns   <- do.call("rbind", by(dat.i, INDICES = dat.i[,group], FUN = function (y) {data.frame ( y[1,group, drop=FALSE], parameter = "Ncases", Mean=length(unique(y[,ID])), stringsAsFactors = FALSE) }))
     tabs <- plyr::rbind.fill(tabs, Ns)
     if(!is.null(group.differences.by))   {
         m    <- tabs
         m$comb.group <- apply(m, 1, FUN = function (ii) { eatTools::crop(paste( ii[group], collapse = "."))})
         m$all.group  <- 1
         tempR<- res.group <- setdiff(group, group.differences.by)
         if(length(res.group) == 0 ) {res.group <- "all.group"}
         difs <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                 if(length(tempR)>0) {
                     datSel <- merge(dat.i, iii[!duplicated(iii[,res.group]),res.group,drop=FALSE], by = res.group, all = FALSE)
                 }  else  {
                     datSel <- dat.i
                 }
                 tbl    <- table(datSel[,c(group.differences.by, dependent)])
                 chisq  <- chisq.test(tbl, correct = correct)
                 scumm  <- iii[!duplicated(iii[,res.group]),res.group,drop = FALSE]
                 group  <- paste( paste( colnames(scumm), as.character(scumm[1,]), sep="="), sep="", collapse = ", ")
                 dif.iii<- data.frame(group = group, parameter = "chiSquareTest", comparison = "groupDiff", depVar = dependent, modus=modus, coefficient = c("chi2","df","pValue"), value = c(chisq[["statistic"]],chisq[["parameter"]],chisq[["p.value"]]) , stringsAsFactors = FALSE )
                 return(dif.iii)}))                                             ### siehe http://www.vassarstats.net/dist2.html
     }                                                                          ### und http://onlinestatbook.com/2/tests_of_means/difference_means.html
     ret  <- reshape2::melt(tabs, measure.vars = c("Mean", "std.err"), na.rm=TRUE)
     ret[,"coefficient"] <- car::recode(ret[,"variable"], "'Mean'='est'; 'std.err'='se'")
     ret  <- data.frame ( group = apply(ret[,group,drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = dependent, modus = modus, comparison = NA, ret[,c("coefficient", "parameter")], value = ret[,"value"], ret[,group,drop=FALSE], stringsAsFactors = FALSE)
     if(!is.null(group.differences.by))   {return(eatTools::facToChar(plyr::rbind.fill(ret,difs)))} else {return(eatTools::facToChar(ret))}}


jackknife.table <- function ( dat.i , a) {
                   for ( i in names(a)) { assign(i, a[[i]]) }
                   dat.i[,dependent] <- factor(dat.i[,dependent], levels = expected.values)
                   typeS     <- car::recode(type, "'JK2'='JKn'")
                   design    <- survey::svrepdesign(data = dat.i[,c(group, dependent)], weights = dat.i[,wgt], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(dat.i[,ID], repA[,ID] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                   formel    <- as.formula(paste("~factor(",dependent,", levels = expected.values)",sep=""))
                   means     <- survey::svyby(formula = formel, by = as.formula(paste("~", paste(as.character(group), collapse = " + "))), design = design, FUN = svymean, deff = FALSE, return.replicates = TRUE)
                   Ns        <- do.call("rbind", by(dat.i, INDICES = dat.i[,group], FUN = function (y) {data.frame ( y[1,group, drop=FALSE], variable = "est____________Ncases", value=length(unique(y[,ID])), stringsAsFactors = FALSE) }))
                   cols      <- match(paste("factor(",dependent,", levels = expected.values)",expected.values,sep=""), colnames(means))
                   colnames(means)[cols] <- paste("est",expected.values, sep="____________")
                   cols.se   <- grep("^se[[:digit:]]{1,5}$", colnames(means) )
                   stopifnot(length(cols) == length(cols.se))
                   colnames(means)[cols.se] <- paste("se____________", expected.values, sep="")
                   molt      <- reshape2::melt(data=means, id.vars=group, na.rm=TRUE)
                   molt      <- rbind(molt, Ns)
                   splits    <- data.frame ( do.call("rbind", strsplit(as.character(molt[,"variable"]),"____________")), stringsAsFactors = FALSE)
                   colnames(splits) <- c("coefficient", "parameter")
                   ret       <- data.frame ( group = apply(molt[,group,drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = dependent, modus = paste(modus,"survey", sep="__"), comparison = NA, splits, value = molt[,"value"], molt[,group,drop=FALSE], stringsAsFactors = FALSE)
                   if(!is.null(group.differences.by))   {
                      m            <- ret
                      m$comb.group <- apply(m, 1, FUN = function (ii) { eatTools::crop(paste( ii[group], collapse = "."))})
                      m$all.group  <- 1
                      res.group    <- tempR <- setdiff(group, group.differences.by)
                      if(length(res.group) == 0 ) {res.group <- "all.group"}
                      difs           <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                        if(length(tempR)>0) {
                                           datSel    <- merge(dat.i, iii[!duplicated(iii[,res.group]),res.group,drop=FALSE], by = res.group, all = FALSE)
                                        }  else  {
                                           datSel    <- dat.i
                                        }
                                        designSel <- survey::svrepdesign(data = datSel[,c(group.differences.by, dependent)], weights = datSel[,wgt], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(datSel[,ID], repA[,ID] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                                        formel    <- as.formula(paste("~", group.differences.by ,"+",dependent,sep=""))
                                        tbl       <- survey::svychisq(formula = formel, design = designSel, statistic = "Chisq")
                                        scumm     <- iii[!duplicated(iii[,res.group]),res.group,drop = FALSE]
                                        group     <- paste( paste( colnames(scumm), as.character(scumm[1,]), sep="="), sep="", collapse = ", ")
                                        dif.iii   <- data.frame(group = group, parameter = "chiSquareTest", depVar = dependent, modus = paste(modus,"survey", sep="__"), coefficient = c("chi2","df","pValue"), value = c(tbl[["statistic"]],tbl[["parameter"]],tbl[["p.value"]]) , stringsAsFactors = FALSE )
                                        return(dif.iii)                         ### siehe http://www.vassarstats.net/dist2.html
                      } ))                                                      ### http://onlinestatbook.com/2/tests_of_means/difference_means.html
                      difs[,"comparison"] <- "groupDiff"
                   }
                   if(!is.null(group.differences.by))   {return(eatTools::facToChar(plyr::rbind.fill(ret,difs)))} else {return(eatTools::facToChar(ret))}}


conv.mean      <- function (dat.i , a) {
                  for ( i in names(a)) { assign(i, a[[i]]) }
                  deskr    <- data.frame ( do.call("rbind", by(data = dat.i, INDICES = dat.i[,group], FUN = function ( sub.dat) {
                              prefix <- sub.dat[1,group, drop=FALSE]
                              if ( all(sub.dat[,wgt] == 1) )  {
                                   useWGT <- NULL
                              }  else  {
                                   useWGT <- sub.dat[,wgt]
                                   attr(useWGT, "onlyUnweightedN") <- TRUE      ### das ist noetig, weil jetzt immer ungewichtete Ns bestimmt werden sollen,
                              }                                                 ### auch wenn der user gewichte angibt. Das wird ja an decr() uebergeben
                              ret    <- data.frame ( nValidUnweighted = length(na.omit(sub.dat[, dependent ])), prefix, eatTools::descr(sub.dat[, dependent ], p.weights = useWGT, na.rm=na.rm)[,c("N", "N.valid", "Mean", "std.err", "Var", "SD")], stringsAsFactors = FALSE)
                              names(ret) <- c( "nValidUnweighted", group , "Ncases", "NcasesValid", "mean", "se.mean", "var","sd")
                              return(ret)})), modus=modus, stringsAsFactors = FALSE)
                  if(!is.null(group.differences.by))   {
                     nCat <- table(as.character(dat.i[,group.differences.by]))
                     if ( length(nCat) < 2 ) {
                          cat(paste("Warning: Grouping variable '", group.differences.by, "' only has one category within imputation and/or nest. Group differences cannot be computed. Skip computation.\n",sep=""))
                     }  else  {
                          m            <- deskr
                          m$comb.group <- apply(m, 1, FUN = function (ii) { eatTools::crop(paste( ii[group], collapse = "."))})
                          m$all.group  <- 1
                          res.group    <- tempR <- setdiff(group, group.differences.by)
                          if(length(res.group) == 0 ) {res.group <- "all.group"}
                          kontraste    <- combinat::combn ( x = sort(unique(as.character(m[,group.differences.by]))), m = 2, simplify = FALSE)
                          difs         <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                          ret <- do.call("rbind", lapply(kontraste, FUN = function ( k ) {
                                                 if ( sum ( k %in% iii[,group.differences.by]) != length(k) ) {
                                                    warning("Cannot compute contrasts for 'group.differences.by = ",group.differences.by,"'.", immediate. = TRUE)
                                                    return(NULL)
                                                 }  else  {
                                                    vgl.iii   <- iii[iii[,group.differences.by] %in% k ,]
                                                    true.diff <- computeTrueDiffAndOtherDiffs(vgl.iii, dat = dat.i, group.differences.by=group.differences.by, kontr = k, value="mean")[["true"]]
                                                    scumm     <- sapply(vgl.iii[,res.group,drop = FALSE], as.character)
                                                    group     <- paste( paste( colnames(scumm), scumm[1,], sep="="), sep="", collapse = ", ")
                                                    dummy     <- do.call("cbind", lapply ( a%$$%group, FUN = function ( gg ) {
                                                                 ret <- data.frame ( paste ( rev(unique(vgl.iii[,gg])), collapse = " - "))
                                                                 colnames(ret) <- gg
                                                                 return(ret)}))
                                                    dif.iii   <- data.frame(dummy, group = paste(group, paste(k, collapse = ".vs."),sep="____"), parameter = "mean", coefficient = c("est","se"), depVar = dependent, modus=modus, value = c(true.diff, sqrt( sum(vgl.iii[,"sd"]^2 / vgl.iii[,"nValidUnweighted"]) )) , stringsAsFactors = FALSE )
                                                    stopifnot(nrow(dif.iii)==2, nrow(vgl.iii) == 2)
                                                    dummy2    <- dif.iii[1,]
                                                    dummy2[,"coefficient"] <- "es"
                                                    dummy2[,"value"]       <- dif.iii[which(dif.iii[,"coefficient"] == "est"),"value"] / sqrt(0.5*sum(vgl.iii[,"sd"]^2))
                                                    return(rbind(dif.iii, dummy2))# siehe http://www.vassarstats.net/dist2.html
                                                 } }))                          ### http://onlinestatbook.com/2/tests_of_means/difference_means.html
                                          return(ret)}))
                        difs[,"comparison"] <- "groupDiff"
                     }
                  }
                  deskrR   <- reshape2::melt(data = deskr, id.vars = group, measure.vars = setdiff(colnames(deskr), c("nValidUnweighted", "modus", group) ), na.rm=TRUE)
                  deskrR[,"coefficient"] <- car::recode(deskrR[,"variable"], "'se.mean'='se';else='est'")
                  deskrR[,"parameter"]   <- gsub("se.mean","mean",deskrR[,"variable"])
                  deskrR   <- data.frame ( group = apply(deskrR[,group,drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = dependent, modus=modus, comparison = NA, deskrR[,c( "parameter", "coefficient", "value", group)], stringsAsFactors = FALSE)
                  if(!is.null(group.differences.by))   {
                      if ( length (nCat ) >1 ) {
                           return(eatTools::facToChar(plyr::rbind.fill(deskrR,difs)))
                      }   else  {
                           return(eatTools::facToChar(deskrR))
                      }
                  }  else {
                      return(eatTools::facToChar(deskrR))
                  } }

computeTrueDiffAndOtherDiffs <- function (difs, repl, dat, kontr, group.differences.by, value) {
          stopifnot ( nrow(difs) == 2 )
          refSeq<- names(table(dat[,group.differences.by]))                     ### immer fokusgruppe MINUS referenzgruppe
          reihe <- match(kontr, refSeq)                                         ### dazu rausfinden, was Referenz ist ...
          trueD <- difs[match(refSeq[max(reihe)],difs[,group.differences.by]),value] - difs[match(refSeq[min(reihe)],difs[,group.differences.by]),value]
          if(!missing(repl)) {otherD<- repl[,refSeq[max(reihe)]] - repl[,refSeq[min(reihe)]]} else {otherD<- NULL}
          return(list(true = trueD, other = otherD))  }

jackknife.adjust.mean <- function (dat.i , a) {
          for ( i in names(a)) { assign(i, a[[i]]) }
          typeS<- car::recode(type, "'JK2'='JKn'")
          repl <- repA[ match(dat.i[,ID], repA[,ID]),]
          des  <- survey::svrepdesign(data = dat.i[,c(group, dependent, adjust)], weights = dat.i[,wgt], type=typeS, scale = 1, rscales = 1, repweights = repl[,-1, drop = FALSE], combined.weights = TRUE, mse = TRUE, rho=rho)
          if ( useEffectLiteR ) {
               ret <- survey::withReplicates(des, funAdjustEL, allNam=a[allNam], return.replicates=TRUE)
          }  else  {
               ret <- survey::withReplicates(des, funAdjust, allNam=a[allNam], return.replicates=TRUE)
          }
          rs   <- m <- data.frame ( group = rep(rownames(as.data.frame ( ret)),2) , depVar = dependent, modus = paste(modus, "survey", sep="__"), comparison = NA, parameter = "mean", coefficient = rep(c("est", "se"), each = nrow(as.data.frame (ret))), value = reshape2::melt(as.data.frame ( ret), measure.vars = colnames(as.data.frame ( ret)))[,"value"], rbind(do.call("rbind",  by(data=dat.i, INDICES = dat.i[,group], FUN = function ( x ) { x[1,group, drop=FALSE]}, simplify = FALSE)),do.call("rbind",  by(data=dat.i, INDICES = dat.i[,group], FUN = function ( x ) { x[1,group, drop=FALSE]}, simplify = FALSE))), stringsAsFactors=FALSE)
          if(!is.null(group.differences.by))   {
             nCat <- table(as.character(dat.i[,group.differences.by]))
             if ( length(nCat) < 2 ) {
                  warning("Grouping variable '", group.differences.by, "' only has one category within imputation and/or nest. Group differences cannot be computed. Skip computation.", immediate. = TRUE)
             }  else  {
                m$comb.group <- apply(m, 1, FUN = function (ii) {eatTools::crop(paste( ii[group], collapse = "."))})
                repl1<- data.frame ( repl = rownames(ret[["replicates"]]),ret[["replicates"]], stringsAsFactors=FALSE)
                m$all.group    <- 1
                res.group      <- tempR  <- setdiff(group, group.differences.by)
                if(length(res.group) == 0 ) {res.group <- "all.group"}
                kontraste      <- combinat::combn ( x = sort(unique(as.character(m[,group.differences.by]))), m = 2, simplify = FALSE)
                difs           <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                  ret <- do.call("rbind", lapply(kontraste, FUN = function ( k ) {
                                         if ( sum ( k %in% iii[,group.differences.by]) != length(k) ) {
                                              warning("Cannot compute contrasts for 'group.differences.by = ",group.differences.by,"'.", immediate. = TRUE)
                                              return(NULL)                      ### Quelle fuer dieses Vorgehen:
                                         } else {                               ### Mail SW an ZKD, 07.11.2012, 17.54 Uhr, "in Absprache mit Dirk"
                                              vgl.iii <- iii[iii[,group.differences.by] %in% k ,]
                                              vgl.iii <- vgl.iii[which(vgl.iii[,"coefficient"] == "est"),]
                                              diffs   <- computeTrueDiffAndOtherDiffs(difs = vgl.iii, repl = repl1, dat=dat.i, kontr = k, group.differences.by=group.differences.by, value="value")
                                              scumm   <- sapply(vgl.iii[,res.group,drop = FALSE], as.character)
                                              group   <- paste( paste( colnames(scumm), scumm[1,], sep="="), sep="", collapse = ", ")
                                              dummy   <- do.call("cbind", lapply ( a%$$%group, FUN = function ( gg ) {
                                                         ret <- data.frame ( paste ( rev(unique(vgl.iii[,gg])), collapse = " - "))
                                                         colnames(ret) <- gg
                                                         return(ret)}))
                                              dif.iii <- data.frame(dummy, group = group, vgl = paste(k, collapse = ".vs."), dif = diffs[["true"]], se =  sqrt(sum((diffs[["true"]] - diffs[["other"]])^2)), stringsAsFactors = FALSE )
                                              return(dif.iii)
                                         } }))
                                  return(ret)}))
                difsL<- data.frame ( comparison = "groupDiff", depVar = dependent, reshape2::melt(data = difs, measure.vars = c("dif", "se") , variable.name = "coefficient" , na.rm=TRUE), modus=paste(modus,"survey", sep="__"), parameter = "mean", stringsAsFactors = FALSE)
                difsL[,"coefficient"] <- car::recode(difsL[,"coefficient"], "'se'='se'; 'es'='es'; else = 'est'")
                difsL[,"group"] <- apply(difsL[,c("group","vgl")],1,FUN = function (z) {paste(z,collapse="____")})
                rs   <- rbind(rs,difsL[,-match("vgl", colnames(difsL))])
             }
          }
          return(eatTools::facToChar(rs)) }

### Hilfsfunktion fuer jackknife.adjust.mean: braucht keine Standardfehler, da die mit jackknife bestimmt werden
funAdjust <- function(w, data, allNam){                                         ### 'w' = weights
       data[,allNam[["wgt"]]] <- w
       frml<- as.formula(paste0(allNam[["dependent"]]," ~ ", paste(allNam[["adjust"]], collapse = " + ")))
       if ( all(data[,allNam[["wgt"]]] == 1) )  {                               ### adjustierte Mittelwerte ungewichtet
            means <- by ( data = data, INDICES = data[,allNam[["group"]]], FUN = function ( gr ) {
                     reg <- lm(frml, data = gr)
                     cof1<- coef(reg)                                           ### untere Zeile: innerhalb der Funktion muss mean(dat...) stehen, nicht mean(gr...), sonst kriegt man nur die unadjustieren Mittelwerte
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
       ret   <- as.vector(means)                                                ### Funktion gibt einen numerischen Vektor zurueck, der Namen hat
       names(ret) <- names(means)
       return(ret)}


### Hilfsfunktion fuer jackknife.adjust.mean
funAdjustEL <- function(w, data, allNam){
       data[,allNam[["wgt"]]] <- w
       if ( all(data[,allNam[["wgt"]]] == 1) )  {                               ### adjustierte Mittelwerte ungewichtet
            res <- EffectLiteR::effectLite(y = allNam[["dependent"]], x = allNam[["group"]], z = allNam[["adjust"]], data = data, fixed.cell = TRUE, fixed.z = FALSE, homoscedasticity = FALSE, method="sem")
       }  else  {
            frml<- as.formula(paste0("~", allNam[["wgt"]]))
            res <- suppressMessages(EffectLiteR::effectLite(y = allNam[["dependent"]], x = allNam[["group"]], z = allNam[["adjust"]], data = data, fixed.cell = TRUE, fixed.z = FALSE, homoscedasticity = FALSE, weights = frml, method="sem"))
       }
       ret <- res@results@adjmeans[,"Estimate"]
       names(ret) <- names(table(data[,allNam[["group"]]]))
       return(ret)}  
       
conv.adjust.mean <- function ( dat.i, a) {
       for ( i in names(a)) { assign(i, a[[i]]) }
       if(isTRUE(useEffectLiteR)) {                                             ### benutze effectLite
           if ( all(dat.i[,wgt] == 1) ) {                                       ### ohne Gewichte
                res <- EffectLiteR::effectLite(y = dependent, x = group, z = adjust, data = dat.i, fixed.cell = TRUE, fixed.z = FALSE, homoscedasticity = FALSE, method="sem")
           }  else  {
                st  <- paste("suppressMessages(EffectLiteR::effectLite(y = dependent, x = group, z = adjust, data = dat.i, fixed.cell = TRUE, fixed.z = FALSE, homoscedasticity = FALSE, weights = ~ ",wgt,", method=\"sem\"))", sep="")
                res <- eval(parse( text=st))
           }
           vals <- reshape2::melt(res@results@adjmeans, measure.vars = c("Estimate", "SE"))[,"value"]
       }  else  {                                                               ### benutze kein effectLite: Standardfehler mit deltamethode
           means <- do.call("rbind", by ( data = dat.i, INDICES = dat.i[,group], FUN = function ( gr ) {
                    frml<- paste0(dependent, " ~ ", paste(adjust, collapse = " + "))
                    if ( all(dat.i[,wgt] == 1) ) {                              ### ohne Gewichte
                         reg <- lm(as.formula(frml), data = gr)                 ### ungewichtete Regression linear
                         x_m <- sapply(dat.i[,adjust], mean)                    ### ungewichtete Mittelwerte der adjustierungsvariablen
                    }  else  {                                                  ### untere Zeile: gewichtete Regression
                         reg <- eval(parse(text = paste0("lm(as.formula(frml), data = gr, weights = ",wgt, ")")))
                         x_m <- unlist(lapply (adjust, FUN = function (u) { Hmisc::wtd.mean(dat.i[,u], weights = dat.i[,wgt])}))
                    }                                                           ### obere Zeile: gewichtete Mittelwerte der adjustierungsvariablen
                    cof1<- coef(reg)
                    adj <- cof1[1] + sum(cof1[-1] * x_m)
                    pars<- c(cof1, x_m)
                    mat <- matrix(0, length(cof1) + length(x_m), length(cof1) + length(x_m))
                    mat[1:length(cof1), 1:length(cof1)] <- vcov(reg)
                    if ( all(dat.i[,wgt] == 1) ) {                              ### ungewichtete Varianz
                         mat[(length(cof1)+1):nrow(mat), (length(cof1)+1):ncol(mat)] <- var(gr[,adjust]) / (nrow(gr)-1)
                    }  else  {                                                  ### gewichtete Varianz
                         mat[(length(cof1)+1):nrow(mat), (length(cof1)+1):ncol(mat)] <- cov.wt(gr[,adjust, drop=FALSE], wt = gr[,wgt], cor = FALSE, center = TRUE)[["cov"]] / (nrow(gr)-1)
                    }
                    frm2<- paste0("~x1 + ", paste(paste("x", 2:length(cof1), sep=""), paste("x", (length(cof1)+1):(length(cof1)+length(x_m)), sep=""), collapse=" + ", sep="*"))
                    se  <- msm::deltamethod(as.formula(frm2), pars, mat)
                    return(data.frame ( mw = adj, se = se, stringsAsFactors = FALSE))}))
           vals  <- reshape2::melt(means, measure.vars = c("mw", "se"))[,"value"]
       }
       rs   <- m <- data.frame ( group = rep(names(table(dat.i[,group])) , 2), depVar = dependent, modus = modus, comparison = NA, parameter = "mean", coefficient = rep(c("est", "se"), each = length(vals)/2), value = vals, rbind(do.call("rbind",  by(data=dat.i, INDICES = dat.i[,group], FUN = function ( x ) { x[1,group, drop=FALSE]}, simplify = FALSE)),do.call("rbind",  by(data=dat.i, INDICES = dat.i[,group], FUN = function ( x ) { x[1,group, drop=FALSE]}, simplify = FALSE))), stringsAsFactors=FALSE)
       if(!is.null(group.differences.by))   {                                   ### jetzt gruppendifferenzen, wenn es welche geben soll
           nCat <- table(as.character(dat.i[,group.differences.by]))
           if ( length(nCat) < 2 ) {
                cat(paste("Warning: Grouping variable '", group.differences.by, "' only has one category within imputation and/or nest. Group differences cannot be computed. Skip computation.\n",sep=""))
           }  else  {
                m$comb.group <- apply(m, 1, FUN = function (ii) { eatTools::crop(paste( ii[group], collapse = "."))})
                m$all.group  <- 1
                res.group    <- tempR <- setdiff(group, group.differences.by)
                if(length(res.group) == 0 ) {res.group <- "all.group"}
                kontraste    <- combinat::combn ( x = sort(unique(as.character(m[,group.differences.by]))), m = 2, simplify = FALSE)
                difs         <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                ret <- do.call("rbind", lapply(kontraste, FUN = function ( k ) {
                                       if ( sum ( k %in% iii[,group.differences.by]) != length(k) ) {
                                            warning("Cannot compute contrasts for 'group.differences.by = ",group.differences.by,"'.", immediate. = TRUE)
                                            return(NULL)
                                       }  else  {
                                            vgl.iii   <- eatTools::makeDataFrame(tidyr::pivot_wider(iii[iii[,group.differences.by] %in% k ,], names_from = "coefficient", values_from = "value"), verbose=FALSE)
                                            true.diff <- computeTrueDiffAndOtherDiffs(vgl.iii, dat = dat.i, group.differences.by=group.differences.by, kontr = k, value="est")[["true"]]
                                            scumm     <- sapply(vgl.iii[,res.group,drop = FALSE], as.character)
                                            group     <- paste( paste( colnames(scumm), scumm[1,], sep="="), sep="", collapse = ", ")
                                            dummy     <- do.call("cbind", lapply ( a%$$%group, FUN = function ( gg ) {
                                                         ret <- data.frame ( paste ( rev(unique(vgl.iii[,gg])), collapse = " - "))
                                                         colnames(ret) <- gg
                                                         return(ret)}))
                                            dif.iii   <- data.frame(dummy, group = paste(group, paste(k, collapse = ".vs."),sep="____"), comparison = "groupDiff", parameter = "mean", coefficient = c("est","se"), depVar = dependent, modus=modus, value = c(true.diff, sqrt(sum(vgl.iii[,"se"]^2))) , stringsAsFactors = FALSE )
                                            stopifnot(nrow(dif.iii)==2, nrow(vgl.iii) == 2)
                                            return(dif.iii)
                                                 } }))
                                return(ret)}))
           }
       }
       if(!is.null(group.differences.by))   {
          if ( length (nCat ) >1 ) {
               return(eatTools::facToChar(plyr::rbind.fill(rs,difs)))
          }   else  {
               return(eatTools::facToChar(rs))
          }
       }  else {
          return(eatTools::facToChar(rs))
       } }


jackknife.mean <- function (dat.i , a) {
          for ( i in names(a)) { assign(i, a[[i]]) }                            ### alles auf den namespace exportieren
          typeS<- car::recode(type, "'JK2'='JKn'")
          repl <- repA[ match(dat.i[,ID], repA[,ID]),]
          des  <- survey::svrepdesign(data = dat.i[,c(group, dependent)], weights = dat.i[,wgt], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repl[,-1, drop = FALSE], combined.weights = TRUE, rho=rho)
          rets <- data.frame ( target = c("Ncases", "NcasesValid", "mean", "var"), FunctionToCall = c(NA,NA,"svymean","svyvar"), formelToCall = c("paste(\"~ \", \"N_weighted\",sep=\"\")","paste(\"~ \", \"N_weightedValid\",sep=\"\")","paste(\"~ \",dependent, sep = \"\")","paste(\"~ \",dependent, sep = \"\")"), naAction = c("FALSE","TRUE","na.rm","na.rm"), stringsAsFactors = FALSE)
          ret  <- apply(rets, 1, FUN = function ( toCall ) {                    ### svyby wird dreimal aufgerufen ...
                  if (is.na(toCall[["FunctionToCall"]])) {                      ### Achtung: N und N.valid sollen jetzt immer ungewichtet bestimmt werden! deshalb wird svy fuer diese ersten beiden nicht mehr gecallt
                      resL <- do.call("rbind", by(data=dat.i, INDICES = dat.i[,group], FUN = function (y){
                              if (toCall[["target"]] == "Ncases") {weg <- 0} else {weg <- length(which(is.na(y[,dependent])))}
                              r1 <- data.frame ( y[1,group, drop=FALSE], parameter =toCall[["target"]], coefficient="est", value=nrow(y)-weg, stringsAsFactors = FALSE)
                              return(r1) }))
                  }  else  {
                      do   <- paste("svyby(formula = as.formula(",toCall[["formelToCall"]],"), by = as.formula(paste(\"~\", paste(group, collapse = \" + \"))), design = des, FUN = ",toCall[["FunctionToCall"]],",na.rm=",toCall[["naAction"]],", deff = FALSE, return.replicates = TRUE)",sep="")
                      res  <- eval(parse(text=do)) |> suppressWarnings()        ### Warning erklaert in Word-Doc, wird unterdrueckt da irrelevant fuer Paket
                      resL <- reshape2::melt( data = res, id.vars = group, variable.name = "coefficient" , na.rm=TRUE)
                      stopifnot(length(table(resL[,"coefficient"])) == 2)
                      resL[,"coefficient"] <- car::recode(resL[,"coefficient"], "'se'='se'; else ='est'")
                      resL[,"parameter"]   <- toCall[["target"]]
                      attr(resL, "original") <- res
                  }
                  return(resL)})
          sds  <- do.call("rbind", by(data = dat.i, INDICES =  dat.i[,group], FUN = function (uu) {
                  namen   <- uu[1, group, drop=FALSE]                           ### Warning erklaert in Word-Doc, wird unterdrueckt da irrelevant fuer Paket
                  sub.rep <- repl[ match(uu[,ID], repl[,ID] ) ,  ]
                  des.uu  <- survey::svrepdesign(data = uu[,c(group, dependent)], weights = uu[,wgt], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = sub.rep[,-1, drop = FALSE], combined.weights = TRUE, rho=rho)
                  var.uu  <- survey::svyvar(x = as.formula(paste("~",dependent,sep="")), design = des.uu, deff = FALSE, return.replicates = TRUE, na.rm = na.rm) |> suppressWarnings()
                  ret     <- data.frame(namen, est = as.numeric(sqrt(coef(var.uu))), se =  as.numeric(sqrt(vcov(var.uu)/(4*coef(var.uu)))), stringsAsFactors = FALSE )
                  return(ret)}) )
          sds  <- data.frame ( reshape2::melt(data = sds, id.vars = group, variable.name = "coefficient" , na.rm=TRUE), parameter = "sd", stringsAsFactors = FALSE)
          resAl<- rbind(do.call("rbind",ret), sds)
          resAl<- data.frame ( group = apply(resAl[,group,drop=FALSE],1,FUN = function (z) {paste(z,collapse=group.delimiter)}), depVar = dependent, modus=paste(modus,"survey", sep="__"), comparison = NA, resAl[,c("parameter","coefficient","value",group)] , stringsAsFactors = FALSE)
          if(!is.null(group.differences.by))   {
             nCat <- table(as.character(dat.i[,group.differences.by]))
             if ( length(nCat) < 2 ) {
                  warning("Grouping variable '", group.differences.by, "' only has one category within imputation and/or nest. Group differences cannot be computed. Skip computation.", immediate. = TRUE)
             }  else  {
                m1   <- attr(ret[[ which(rets[,"target"] == "mean") ]], "original")
                m1   <- merge(m1, unique(resAl[,intersect(colnames(resAl), colnames(m1)), drop=FALSE]), all.x = TRUE, all.y = FALSE)
                sd   <- sds[which(sds[,"coefficient"] == "est"),]
                colnames(sd) <- car::recode(colnames(sd), "'value'='standardabweichung'")
                m    <- merge(m1, sd, by = group, all = TRUE)
                stopifnot(nrow(m) == nrow(m1))
                m$comb.group <- apply(m, 1, FUN = function (ii) {eatTools::crop(paste( ii[group], collapse = "."))})
                repl1<- data.frame(t(attr(attr(ret[[which(rets[,"target"] == "mean")]], "original"), "replicates")), stringsAsFactors = FALSE )
                colnames(repl1) <- replCols <- paste("replNum", 1:ncol(repl1), sep="")
                repl1[,"comb.group"] <- rownames(repl1)
                m    <- merge(m, repl1, by = "comb.group" )
                m$all.group    <- 1
                res.group      <- tempR  <- setdiff(group, group.differences.by)
                if(length(res.group) == 0 ) {res.group <- "all.group"}
                kontraste      <- combinat::combn ( x = sort(unique(as.character(m[,group.differences.by]))), m = 2, simplify = FALSE)
                difs           <- do.call("rbind", by(data = m, INDICES = m[,res.group], FUN = function (iii)   {
                                  ret <- do.call("rbind", lapply(kontraste, FUN = function ( k ) {
                                         if ( sum ( k %in% iii[,group.differences.by]) != length(k) ) {
                                              warning("Cannot compute contrasts for 'group.differences.by = ",group.differences.by,"'.", immediate. = TRUE)
                                              return(NULL)                      ### Quelle fuer dieses Vorgehen:
                                         } else {                               ### Mail SW an ZKD, 07.11.2012, 17.54 Uhr, "in Absprache mit Dirk"
                                              vgl.iii <- iii[iii[,group.differences.by] %in% k ,]
                                              repl    <- eatTools::makeDataFrame(t(vgl.iii[,replCols]), verbose=FALSE)
                                              colnames(repl) <- vgl.iii[,group.differences.by]
                                              diffs   <- computeTrueDiffAndOtherDiffs(difs = vgl.iii, repl =  repl, dat=dat.i, kontr = k, group.differences.by=group.differences.by, value=intersect(colnames(vgl.iii), dependent))
                                              scumm   <- sapply(vgl.iii[,res.group,drop = FALSE], as.character)
                                              group   <- paste( paste( colnames(scumm), scumm[1,], sep="="), sep="", collapse = ", ")
                                              dummy   <- do.call("cbind", lapply ( a%$$%group, FUN = function ( gg ) {
                                                         ret <- data.frame ( paste ( rev(unique(vgl.iii[,gg])), collapse = " - "))
                                                         colnames(ret) <- gg    ### obere Zeile: Achtung! Hier muss a%$$%group stehen, nicht group!!
                                                         return(ret)}))
                                              dif.iii <- data.frame(dummy, group = group, vgl = paste(k, collapse = ".vs."), dif = diffs[["true"]], se =  sqrt(sum((diffs[["true"]] - diffs[["other"]])^2)), stringsAsFactors = FALSE )
                                              dif.iii <- data.frame(dif.iii, es = dif.iii[["dif"]] / sqrt(0.5*sum(vgl.iii[,"standardabweichung"]^2)))
                                              return(dif.iii)
                                         } }))
                                  return(ret)}))
                difsL<- data.frame ( depVar = dependent, reshape2::melt(data = difs, measure.vars = c("dif", "se", "es") , variable.name = "coefficient" , na.rm=TRUE), modus=paste(modus,"survey", sep="__"), parameter = "mean", stringsAsFactors = FALSE)
                difsL[,"coefficient"] <- car::recode(difsL[,"coefficient"], "'se'='se'; 'es'='es'; else = 'est'")
                difsL[,"comparison"]  <- "groupDiff"
                difsL[,"depVar"]      <- dependent
                difsL[,"group"] <- apply(difsL[,c("group","vgl")],1,FUN = function (z) {paste(z,collapse="____")})
                resAl<- rbind(resAl,difsL[,-match("vgl", colnames(difsL))])
             }
          }
          return(eatTools::facToChar(resAl)) }

### richtige Reihenfolge des gepasteten strings im Rueckgabeobjekt
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

jackknife.cov <- function (dat.i , a){
          for ( i in names(a)) { assign(i, a[[i]]) }
          typeS<- car::recode(type, "'JK2'='JKn'")
          repl <- repA[ match(dat.i[,ID], repA[,ID]),]
          des  <- survey::svrepdesign(data = dat.i[,c(group, dependent)], weights = dat.i[,wgt], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repl[,-1, drop = FALSE], combined.weights = TRUE, rho=rho)
          ret  <- survey::withReplicates(des, groupVersusTotalMean, allNam=a[allNam])
          rs   <- data.frame ( group =  buildString(dat= dat.i,allNam=a[allNam], refGrp=refGrp, reihenfolge) , depVar = dependent, modus = NA, comparison = "crossDiff", parameter = "mean", coefficient = rep(c("est", "se"), each = nrow(ret)), value = reshape2::melt(as.data.frame ( ret), measure.vars = colnames(as.data.frame ( ret)))[,"value"], rbind(do.call("rbind",  by(data=dat.i, INDICES = dat.i[,group], FUN = function ( x ) { x[1,group, drop=FALSE]}, simplify = FALSE)),do.call("rbind",  by(data=dat.i, INDICES = dat.i[,group], FUN = function ( x ) { x[1,group, drop=FALSE]}, simplify = FALSE))), stringsAsFactors=FALSE)
          return(rs)}


### Hilfsfunktion fuer jackknife.cov und conv.cov
groupVersusTotalMean <- function(w, data, allNam){                              ### 'd' = dependent; 'g' = grouping; 'w' = weights
       data[,allNam[["wgt"]]] <- w
       if ( all(data[,allNam[["wgt"]]] == 1) )  {
           ret <- by(data=data, INDICES = data[,allNam[["group"]]], FUN = function ( y ) { mean(y[,allNam[["dependent"]]])})
           gm  <- mean(data[,allNam[["dependent"]]])                            ### Gesamtmittelwert ungewichtet
       }  else  {
           ret <- by(data=data, INDICES = data[,allNam[["group"]]], FUN = function ( y ) { Hmisc::wtd.mean(y[,allNam[["dependent"]]], weights = y[,allNam[["wgt"]]])})
           gm  <- Hmisc::wtd.mean(data[,allNam[["dependent"]]], weights = data[,allNam[["wgt"]]])
       }                                                                        ### obere zeile: Gesamtmittelwert gewichtet
       ret <- gm - ret
       return(ret)}

conv.cov <- function (dat.i, a){
          covs<- boot::boot(data=dat.i, R = a%$$%nBoot, statistic = function ( x, i) {groupVersusTotalMean(w = x[i,a%$$%wgt], data = x[i,c(a%$$%group, a%$$%dependent)], allNam=a[a%$$%allNam])})
          mns <- colMeans(covs$t)
          ses <- sapply(as.data.frame(covs$t), FUN = sd)                        
          rs  <- data.frame ( group =  buildString(dat= dat.i,allNam=a[a%$$%allNam], refGrp=a%$$%refGrp, a%$$%reihenfolge) , depVar = a%$$%dependent, modus = NA, comparison = "crossDiff", parameter = "mean", coefficient = rep(c("est", "se"), each = length(mns)), value = c(mns, ses), rbind(do.call("rbind",  by(data=dat.i, INDICES = dat.i[,a%$$%group], FUN = function ( x ) { x[1,a%$$%group, drop=FALSE]}, simplify = FALSE)),do.call("rbind",  by(data=dat.i, INDICES = dat.i[,a%$$%group], FUN = function ( x ) { x[1,a%$$%group, drop=FALSE]}, simplify = FALSE))), stringsAsFactors=FALSE)
          return(rs)}
          
### Hilfsfunktion fuer repGlm()
jackknife.glm <- function (dat.i , a) {
                 for ( i in names(a)) { assign(i, a[[i]]) }                     ### family und formula muessen sowieso auf den NAMESPACE gegeben werden, sonst misslingt svyglm. Keine Ahnung warum.
                 sub.ana <- by(data = dat.i, INDICES = dat.i[,group], FUN = function (sub.dat) {
                            nam    <- sub.dat[1,group,drop=FALSE]               ### Funktionen, die die Liste 'a' nicht mehr weitergeben muessen, duerfen auf den namespace exportieren
                            if ( wgt == "wgtOne") {
                                 glm.ii <- test <- glm(formula = formula, data = sub.dat, family = family)
                            }  else  {
                                 glm.ii <- test <- eval(parse(text = paste("glm(formula = formula, data = sub.dat, family = family, weights = ",wgt,")",sep=""))) |> suppressWarnings()
                            }
                            singular       <- names(glm.ii[["coefficients"]])[which(is.na(glm.ii[["coefficients"]]))]
                            if(doJK) {
                                modus    <- paste(modus, "survey", sep="__")
                                typeS    <- car::recode(type, "'JK2'='JKn'")
                                design   <- survey::svrepdesign(data = sub.dat[,c(group, independent, dependent) ], weights = sub.dat[,wgt], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(sub.dat[,ID], repA[,ID] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                                if(length(singular) == 0 & isFALSE(a%$$%forceSingularityTreatment) ) {
    ### hier koennen die Warnungen unterdrueckt werden, denn sie werden ja bereits oben ausgegeben (hier nur fuer jedes replicate separat)
                                   glm.ii  <- survey::svyglm(formula = formula, design = design, return.replicates = FALSE, family = family) |> suppressWarnings()
                                }
                            }
                            summaryGlm     <- summary(glm.ii)
    ### AIC, deviance etc. nur auslesen, wenn linkfunktion NICHT 'identity'. Achtung das na.omit(glm.ii[["coefficients"]]) macht man, um die singulaeren Terme loszuwerden
                            if (  length(grep("gaussian", glm.ii[["family"]][["family"]])) > 0 ) {
                                  r2     <- list ( R2 = var(glm.ii$fitted.values)/var(glm.ii$y) , N = length(glm.ii$fitted.values) )
                                  res.bl <- data.frame ( group=paste(sub.dat[1,group], collapse=group.delimiter), depVar =dependent,modus = modus, parameter = c(rep(c("Ncases",names(na.omit(glm.ii[["coefficients"]]))),2),"R2"),
                                            coefficient = c(rep(c("est","se"),each=1+length(names(na.omit(glm.ii[["coefficients"]])))),"est"),
                                            value=c(r2[["N"]],na.omit(glm.ii[["coefficients"]]),NA,summaryGlm$coef[,2],r2[["R2"]]),sub.dat[1,group, drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
                            }   else  {                                         ### hier jetzt fuer logistische modelle
                                  r2     <- fmsb::NagelkerkeR2(glm.ii)          ### pseudo-R2 fuer logistische Modelle, herkoemmliches fuer lineare
                                  if(wgt != "wgtOne" && doJK == FALSE) {        ### bug (?) in nagelkerke funktion wenn r2 fuer gewichtete logistische regressionen ohne jackknife bestimmt wird. deswegen wird hier r2 immer ungewichtet bestimmt
                                     unwgt <- glm(formula = formula, data = sub.dat, family = family)
                                     r2    <- fmsb::NagelkerkeR2(unwgt)
                                  }
                                  res.bl <- data.frame ( group=paste(sub.dat[1,group], collapse=group.delimiter), depVar =dependent,modus = modus, parameter = c(rep(c("Ncases",names(na.omit(glm.ii[["coefficients"]]))),2),"R2","deviance", "null.deviance", "AIC", "df.residual", "df.null"),
                                            coefficient = c(rep(c("est","se"),each=1+length(names(na.omit(glm.ii[["coefficients"]])))),rep("est", 6)), value=c(r2[["N"]],na.omit(glm.ii[["coefficients"]]),NA,summaryGlm$coef[,2],r2[["R2"]],test$deviance, test$null.deviance, test$aic, test$df.residual, test$df.null),sub.dat[1,group, drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
                            }
                            attr(glm.ii, "r2") <- r2                            ### r2 ranhaengen, damit es bei reconstructResultsStructureGlm nicht neu berechnet werden muss
                            if(doJK || isTRUE(useWec) ) {                       ### jetzt kommt die Behandlung, wenn zusaetzlich zu JK noch singularitaet auftritt. das ueberschreibt nun die bisherige "res.bl"
                                if(length(which(is.na(glm.ii[["coefficients"]]))) > 0 ) {
                                   message("Singularity problem in regression estimation for ", length(singular)," coefficient(s): ",paste(singular, collapse = ", "),". Try workaround ... ")
                                }
                                if(isTRUE(forceSingularityTreatment) && isFALSE(useWec)) {
                                   message("Compute coefficients in the expectation of singularities ... ")
                                }
                                if(length(singular) > 0 || forceSingularityTreatment == TRUE ) {
                                   stopifnot(length(as.character(formula)) == 3 )
                                   if ( isFALSE(useWec) ) {
                                       warning("Unidentified bug with Nagelkerkes r^2 in singularity treatment. No r^2 is computed.", immediate. = TRUE)
                                       if ( glmTransformation == "none" )  {resRoh <- survey::withReplicates(design, getOutputIfSingular, allNam=a[allNam],  frml= formula, fam =family)}
                                       if ( glmTransformation == "sdY" )   {resRoh <- survey::withReplicates(design, getOutputIfSingularT1, allNam=a[allNam],  frml= formula, fam =family)}
                                   }  else  {                                   ### untere Zeile: Kontraste fuer jedes replicate separat bestimmen
                                       if ( crossDiffSE.engine == "lavaan") {
                                            if(doJK ) {                         ### mit lavaan und mit replikationsanalyse
                                                design1<- survey::svrepdesign(data = dat.i[,c(independent, dependent)], weights = dat.i[,wgt], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(dat.i[,ID], repA[,ID] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                                                ret    <- survey::withReplicates(design1, funadjustLavaanWec, allNam=a[allNam])
                                                resRoh <- data.frame ( ret, stringsAsFactors=FALSE)
                                                rownames(resRoh) <- paste0(independent, rownames(resRoh))
                                            }  else  {                          ### mit lavaan, aber ohne ohne Replikationsanalysen
                                                res    <- groupToTotalMeanComparisonLavaan ( d = dat.i, y.var = dependent, group.var = independent, weight.var=wgt, heterogeneous=hetero, stchgrsz=stochasticGroupSizes, extended.results=FALSE, lavaan.syntax.path=NULL, run=TRUE, lavaan.summary.output=FALSE )
                                                resRoh <- data.frame ( theta = res[,"est"], SE = res[,"se"], stringsAsFactors=FALSE)
                                                rownames(resRoh) <- paste0(independent, res[,"par"])
                                            }
                                       }  else  {                               ### bisherge lm-methode
                                            contrasts(dat.i[,as.character(formula)[3]]) <- eatTools::contr.wec.weighted(dat.i[,as.character(formula)[3]], omitted=names(table(dat.i[,as.character(formula)[3]]))[1], weights = dat.i[,wgt])
                                            if(doJK ) {
                                                design1<- survey::svrepdesign(data = dat.i[,c(independent, dependent) ], weights = dat.i[,wgt], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(dat.i[,ID], repA[,ID] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
    ### fuer Replikationsanalysen werden die Standardfehler nicht gebraucht, sondern ueber jackknife bestimmt. Man braucht nur die Koeffizienten. Das heisst, hierfuer kann die konventionelle lm() Funktion benutzt werden
    ### Achtung: in untere Zeile MUSS data.frame(...) vorgeschrieben werden, weil das Objekt sonst Klasse 'svrepstat' hat und nicht ge-rbinded werden kann
                                                resRoh1<- data.frame ( survey::withReplicates(design1, getOutputIfSingularWec, allNam=a[allNam], frml = formula), stringsAsFactors=FALSE)
                                                contrasts(dat.i[,as.character(formula)[3]]) <- eatTools::contr.wec.weighted(dat.i[,as.character(formula)[3]], omitted=names(table(dat.i[,as.character(formula)[3]]))[length(names(table(dat.i[,as.character(formula)[3]])))], weights =  dat.i[,wgt])
                                                design2<- survey::svrepdesign(data = dat.i[,c(independent, dependent) ], weights = dat.i[,wgt], type=typeS, scale = scale, rscales = rscales, mse=mse, repweights = repA[match(dat.i[,ID], repA[,ID] ),-1,drop = FALSE], combined.weights = TRUE, rho=rho)
                                                resRoh2<- data.frame ( survey::withReplicates(design2, getOutputIfSingularWec, allNam=a[allNam], frml = formula), stringsAsFactors=FALSE)
                                                resRoh <- rbind(resRoh1, resRoh2)[which(!duplicated(c(row.names(resRoh1), row.names(resRoh2)))),]
                                            }  else  {                          ### linear models with weighted observations: https://www.r-bloggers.com/2015/09/linear-models-with-weighted-observations/
                                                res1   <- eval(parse(text=createCall ( hetero=hetero, allNam=a[allNam], formula=formula) ), envir=NULL)
                                                contrasts(dat.i[,as.character(formula)[3]]) <- eatTools::contr.wec.weighted(dat.i[,as.character(formula)[3]], omitted=names(table(dat.i[,as.character(formula)[3]]))[length(names(table(dat.i[,as.character(formula)[3]])))], weights =  dat.i[,wgt])
                                                res2   <- eval(parse(text=createCall ( hetero=hetero, allNam=a[allNam], formula=formula) ))
                                                rowNams<- c(rownames(summary(res1)[["coefficients"]]), rownames(summary(res2)[["coefficients"]]))
                                                weg    <- which(duplicated(rowNams))
                                                rowNams<- rowNams[-weg]         ### Achtung!! Das Problem ist, dass diese zeilen sowohl fuer lm() als auch fuer estimatr::lm_robust funktionieren muessen!
                                                resRoh <- data.frame ( rbind ( summary(res1)[["coefficients"]], summary(res2)[["coefficients"]]), stringsAsFactors = FALSE)[-weg,1:2]
                                                colnames(resRoh) <- c("theta", "SE")
                                                resRoh <- rbind(resRoh, data.frame ( theta = c(summary(res1)[["r.squared"]], length(res1[["fitted.values"]])), SE = c(NA, NA), stringsAsFactors = FALSE))
                                                rownames(resRoh) <- c(rowNams, "R2", "Nvalid")
                                            }
                                        }
                                   }                                            ### rownames(resRoh) <- gsub("N", "Ncases", rownames(resRoh))
                                   index      <- which(nchar(rownames(resRoh)) == 0)
                                   if(length(index)>0) { for ( j in 1:length(index)) { rownames(resRoh)[index[j]] <- paste("dummyPar",j,sep="")}}
                                   res.bl     <- data.frame ( group=paste(sub.dat[1,group], collapse=group.delimiter), depVar =dependent,modus = modus, parameter = rep(rownames(resRoh), 2), coefficient = c(rep("est", nrow(resRoh)), rep("se", nrow(resRoh))),
                                                 value = c(resRoh[,"theta"], resRoh[,"SE"]), sub.dat[1,group, drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
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
                 
### Hilfsfunktion fuer jackknife.glm (wec mit lavaan)
### stochastische gruppengroessen werden nicht parametrisiert, da die ja immer irgendwie nicht stochastisch sind wenn jackknife gemacht wird
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
    ### ggf. nach Trendgruppen getrennt
            trend <- lapply(names(jk2.out[["resT"]]), FUN = function (tr) {
                     cat(paste("       Trend group: '", tr, "'.\n",sep=""))
                     splitData <- by ( data = jk2.out[["resT"]][[tr]], INDICES = jk2.out[["resT"]][[tr]][,c("group", "depVar")], FUN = function ( spl ) {return(spl)})
			               if(is.null(analyses)) {analyses <- 1:length(splitData)}
                     for ( i in analyses) {
                           spl    <- splitData[[i]]
                           weg1   <- which ( spl[,"parameter"] %in% c("Ncases","Nvalid","R2", "deviance", "df.null", "df.residual", "null.deviance", "AIC"))
                           weg2   <- grep("wholePopDiff",spl[,"parameter"])
                           weg3   <- grep("^p$", spl[,"coefficient"])
                           ret    <- reshape2::dcast(spl[-unique(c(weg1, weg2, weg3)),], parameter~coefficient)
                           ret[,"t.value"] <- ret[,"est"] / ret[,"se"]
                           df     <- spl[ spl[,"parameter"] == "Ncases" & spl[,"coefficient"] == "est"  ,"value"] - nrow(ret)
                           ret[,"p.value"] <- 2*(1-pt( q = abs(ret[,"t.value"]), df = df ))
                           ret[,"sig"]     <- eatTools::num.to.cat(x = ret[,"p.value"], cut.points = c(0.001, 0.01, 0.05, 0.1), cat.values = c("***", "**", "*", ".", ""))
                           retNR  <- ret
                           ret    <- eatTools::roundDF(ret, digits=digits)
                           groupNamen <- setdiff(colnames(spl), c("group","depVar","modus", "parameter", "coefficient","value", "comparison", "id", "row", "type"))
                           if ( length(groupNamen)>0) {
                                cat ( paste( "            groups: ", paste( groupNamen, unlist(lapply(spl[1,groupNamen], as.character)), sep=" = ", collapse = "; "),"\n",sep=""))
                           }
                           if ( length(add) >0 ) {
                                for ( jj in names(add)) {
                                     lz  <- 18 - nchar(jj)                      ### Anzahl leerzeichen
                                     if ( lz < 0 ) {lz <- 0}
                                     cat(paste( paste(rep(" ", times = lz),collapse=""), jj, ": ", add[[jj]], "\n", sep=""))
                                }
                           }
                           cat ( paste( "dependent Variable: '", as.character(spl[1,"depVar"]), "'\n \n", sep=""))
                           print(ret)
                           r2     <- spl[ spl[,"parameter"] == "R2" ,"value"]
                           cat(paste("\n            R-squared: ",round(r2[1],digits = digits),"; SE(R-squared): ",round(r2[2],digits = digits),"\n",sep=""))
                           if ( isTRUE(printDeviance) ) {
                                for ( prms in c("deviance", "null.deviance", "df.null", "df.residual", "AIC")) {
                                      if ( prms %in% spl[,"parameter"] ) {
                                           cat(paste  ( paste(rep(" ", times = 21 - nchar(prms)),collapse=""),prms,": ", round ( spl[intersect(which(spl[,"parameter"] == prms), which(spl[,"coefficient"] == "est")),"value"],digits = digits), "\n", sep=""))
                                      }
                                }
                           }
                           nn     <- spl[ spl[,"parameter"] == "Ncases" & spl[,"coefficient"] =="est" ,"value"]
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

### Hilfsfunktion fuer repGlm() wenn Regression singulaere Terme enthaelt
### Achtung: negelkerke wird irgendwie falsch berechnet, keine Ahnung woran es liegt, wird ausgeschlossen
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
                       
### wie oben, nur werden hier lineare Transformationen der Regressionskoeffizienten erlaubt
### unstandardisierte Koeffizienten werden genutzt, um den logit vorherzusagen
### jede person hat auf dieser logitvariablen einen wert
### standardabweichung der vorhergesagten logits + (pi^2) / 3 wird genutzt, um varianz der latenten variable zu bestimmen
### unstandardisierten logits werden durch varianz der lat. variablen geteilt, um standardisierte zu erhalten
### Achtung: negelkerke wird irgendwie falsch berechnet, keine Ahnung woran es liegt, wird ausgeschlossen
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
                 stopifnot(length(unique(jk2[,"depVar"]))==1)                   ### jetzt die Ncases dazu
                 Nc  <- repTable.output[intersect(which(repTable.output[,"parameter"] %in% c("Ncases", "NcasesValid")),  which(repTable.output[,"coefficient"] == "est")),]
                 if(!is.null(depVarOri)) {                                      ### obere Zeile subset(...) geht nicht in Paketen, siehe Mail Benjamin, 14.11.2022, 14.28 Uhr
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

### definiert die Funktion (mean, table, glm ... ) und die Methode (JK2, BRR, oder konventionell) fuer Funktion 'eatRep'
chooseFunction <- function (datI, a, pb) {
        pb$tick(); flush.console()
        glms  <- grps <- nams <- NULL                                           ### glms initialisieren, wird ueberschrieben, wenn toCall = 'glm'
        if( a%$$%toCall != "glm" ) {
            fun   <- gsub("\\.\\.", ".", paste(car::recode(a%$$%doJK, "1='jackknife'; 0='conv'"), ifelse(is.null(a%$$%adjust), "", "adjust"), a%$$%toCall, sep="."))
            ana.i <- eval(parse(text=paste0(fun, "(dat.i = datI, a=a)")))
        }  else  {                                                              ### initiate checks specifically for regression models
            doChek<- checkRegression ( dat = datI, allNam=a[a%$$%allNam], useWec=a%$$%useWec)
            ana.i <- jackknife.glm ( dat.i = datI , a=a)
            glms  <- ana.i[["ori.glm"]]
            grps  <- ana.i[["nams"]]
            nams  <- lapply(grps, FUN = function ( x ) { paste(as.character(unlist(x)), collapse=a%$$%group.delimiter)})
            ana.i <- ana.i[["sub.ana"]]
        }
        ana.i <- data.frame ( ana.i, datI[1,c(a%$$%nest, a%$$%imp),drop=FALSE], stringsAsFactors = FALSE, row.names = NULL)
        return(list(ana.i=ana.i, glms=glms, nams=nams, grps=grps))}


findEstSeNames <- function (out) {
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
        return(c(est, se))}

reconstructResultsStructureGlm <- function ( group, neu, grps, group.delimiter, pooled, allNam, modus, formula) {
        r2  <- lapply(neu[[group]], FUN = function ( x ) {attr(x, "r2")[["R2"]]})
        Nval<- lapply(neu[[group]], FUN = function ( x ) {length(x$fitted.values)})
        mat <- which ( unlist(lapply(grps[[1]], FUN = function ( x ) { paste(as.character(unlist(x)), collapse = group.delimiter)})) == group)
    ### achtung: manchmal heissen die Spalten im glm-Output 'est', manchmal 'estimate' ... bzw. manchmal 'se', manchmal 'std.error'
        out <- summary(pooled[[group]])
        nams<- findEstSeNames(out)
    ### in alten Versionen von mice/miceadds waren die Parameter die row.names des 'out'-Objekts
    ### inzwischen sind sie in separaten spalten ... rausfinden, in welcher Weise die Parameter enthalten sind
    ### zur sicherheit wird versucht, die Parameternamen aus dem formula-Objekt zu extrahieren
        prms<- all.vars(formula)[-1]
        col <- which(unlist(lapply(out, FUN = function ( co ) { length(unlist(lapply(prms, FUN = function ( l ) { grep(l, co)})))}))>0)
        stopifnot(length(col) <= 1)
        if ( length(col) == 1) { prm <- as.character(out[,col])}
        if ( length(col) == 0) { prm <- rownames(out)}
        ret <- data.frame ( group=group, depVar =allNam[["dependent"]],modus = modus, comparison = NA, parameter = c(rep(c("Ncases","Nvalid",prm),2),"R2"),
               coefficient = c(rep(c("est","se"),each=2+length(prm)),"est") , value= c(NA, min(unlist(Nval)),out[,nams[1]], NA, NA, out[,nams[2]], pool.R2(unlist(r2), unlist(Nval), verbose = FALSE )[["m.pooled"]]), grps[[1]][[mat]], stringsAsFactors = FALSE, row.names = NULL)
        return(ret)}

### MH wollte in grauer Vorzeit, dass unten stehendes keinen Fehler ergibt; da das aber spaeter schiefgeht, muss es hier eine Fehlermeldung geben
### 06.12.2019, SW: habe jetzt wieder Fehlermeldungen draus gemacht
### 27.08.2024, SW. Funktion hat keine Rueckgabe (mehr), checkt nur
checkData <- function ( sub.dat, a) {
        for ( i in names(a)) { assign(i, a[[i]]) }
        if(!is.null(PSU)) {
            nJkZones <- length(table(as.character(sub.dat[,PSU])))
            if(nJkZones<2)  {
               warning("Found group(s) with less than 2 PSUs. Please check your data!", immediate. = TRUE)
            }
        }                                                                       ### untere Zeile: prueft; es darf GAR KEINE Missings geben
        if( (toCall == "table" & isFALSE(separate.missing.indicator)) | (toCall %in% c("mean", "quantile", "glm") & isFALSE(na.rm) ) )  {
            nObserved <- length(which(is.na(sub.dat[, dependent])))
            if(nObserved>0) {                                                   
               if ( toCall %in% c("mean", "quantile", "glm") ) {
                    stop(paste0("Found unexpected missing data in dependent variable '",dependent,"' for at least one group. Please check your data or set 'na.rm==TRUE'.\n"))
               }  else  {
                    warning(paste0("Found unexpected missing data in dependent variable '",dependent,"' for at least one group although 'separate.missing.indicator' was set to 'FALSE'."), immediate. = TRUE)
               }
            }
        }
        nMissing <- length(which(is.na(sub.dat[, dependent])))                  ### es darf NICHT ALLES missing sein
        if(nMissing == nrow(sub.dat))  { stop(paste0("(At least) some groups without any observed data in dependent variable '",dependent,"'. Please check your data!\n"))}}

checkNests <- function (x, allNam, toAppl, gr) {
        if(length(x[,allNam[["ID"]]]) != length(unique(x[,allNam[["ID"]]])))  {stop("ID variable '",allNam[["ID"]],"' is not unique within nests and imputations.")}
        if( length(toAppl[[gr]])>0) { ret <- lapply( toAppl[[gr]], FUN = function ( y ) {table(x[,y])}) } else {ret <- 1}
        return(list(ret=ret)) }


doSurveyAnalyses <- function (datL1, a) {
        nrep<- table(datL1[, c(a%$$%nest, a%$$%imp)])
        nrep<- prod(dim(nrep))
    ### progress bar fuer analysen anzeigen
        cri1<- nrep > 4 & length(unique(datL1[,a%$$%ID]))>2000 & a%$$%doJK
        cri2<- FALSE
        if ( isFALSE(a%$$%doJK)) {
             cri2<- nrep > 9 & length(unique(datL1[,a%$$%ID]))>5000
        }
        if ( a%$$%progress && (isTRUE(cri1) | isTRUE(cri2)) ) {
             pb  <- progress::progress_bar$new( format = paste0(a%$$%str1, " [:bar] :percent in :elapsed"), incomplete = " ", total = nrep, clear = FALSE, width= 75, show_after = 0.01)
        }  else  {
             pb <- list()
             pb$tick <- function (){return(NULL)}
        }
        ana <- do.call("rbind", by(data = datL1, INDICES = datL1[,a%$$%nest], FUN = function ( datN ) {
    ### optional: multicore handling, um die einzelnen Imputationen parallel zu verarbeiten
               if ( !is.null(a%$$%nCores) && a%$$%nCores > 1 && !is.null(a%$$%imp) && length(unique(datN[,a%$$%imp])) > 1) {
                    if ( length(unique(datN[,a%$$%imp])) < a%$$%nCores) {a$nCores <- length(unique(datN[,a%$$%imp]))}
               }  else {
                    a["nCores"] <- list(NULL)
               }                                                                ### Aufruf wird nicht mehr mit allen Argumenten hintereinander copy/paste eingefuegt
               if (is.null(a%$$%nCores)) {                                      ### hier beginnt der single core Fall
                    anaI <- by(data = datN, INDICES = datN[,a%$$%imp], FUN = chooseFunction, a=a, pb=pb)
               }  else  {                                                       ### jetzt multicore
                    doIt <- function (laufnummer,  a, pb ) {
                            if(!exists("repMean"))  { library(eatRep) }
                            dat <- datN[which(datN[,a%$$%imp] == laufnummer),]
                            ret <- chooseFunction ( datI = dat, a=a, pb=pb)
                            return(ret)}
                    beg  <- Sys.time()
                    cl   <- future::makeClusterPSOCK(a%$$%nCores, verbose=FALSE)
                    anaI <- parallel::clusterApply(cl = cl, x = 1:length(unique(datN[,a%$$%imp])), fun = doIt, a=a, pb=pb)
                    names(anaI) <- unique(datN[,a%$$%imp])
                    stopCluster(cl)
                    tme  <- eatTools::removePattern(capture.output(print(Sys.time() - beg, digits=3)), pattern="Time difference of ")
                    message(paste0("Multicore processing of '",a%$$%modus,"', using ",length(unique(datN[,a%$$%imp]))," imputations and ",a%$$%nCores," cores: ",tme))
               }
               glms <- lapply(anaI, FUN = function ( x ) { x[["glms"]]})
               nams <- lapply(anaI, FUN = function ( x ) { x[["nams"]]})
               grps <- lapply(anaI, FUN = function ( x ) { x[["grps"]]})
    ### hier werden jetzt, sofern nicht genestet imputiert wurde, Regressionskoeffizienten nach Methode "mice" gepoolt
    ### die Listenstruktur des Objekts muss umgedreht werden
               if ( a%$$%toCall == "glm" && a%$$%poolMethod == "mice" && length(glms) > 1) {
                    aussen <- unlist(nams[[1]])
                    innen  <- names(nams)
                    for ( j in 1:length(glms)) { names(glms[[j]]) <- aussen }
                    neu    <- lapply(aussen, FUN = function ( a ) {lapply(innen, FUN = function ( i ) { glms[[i]][[a]]})})
                    pooled <- lapply(neu, FUN = function ( x ) { mice::pool(mice::as.mira(x)) } )
                    names(pooled) <- names(neu) <- aussen
    ### jetzt muss hier die ergebnisstruktur rekonstruiert werden und 'anaI' genannt werden
                    anaI   <- do.call("rbind", lapply(aussen, FUN = reconstructResultsStructureGlm, neu=neu, grps=grps, group.delimiter=a%$$%group.delimiter, pooled=pooled, allNam=a[a%$$%allNam], modus=a%$$%modus, formula=a%$$%formula))
               }  else  {                                                       ### bei methode 'scalar' wird erst spaeter gepoolt
                    anaI <- do.call("rbind", lapply(anaI, FUN = function ( x ) { x[["ana.i"]]}))
               }
               return(anaI)}))
        if ( length(unique(ana[,"modus"])) >1 ) {
             warning("Heterogeneous mode: '", paste(unique(ana[,"modus"]), collapse="', '"),"'", immediate. = TRUE)
        }
        mod <- unique(ana[,"modus"])
    ### es wird nur gepoolt, wenn es mehr als eine Imputation gibt!
        if ( a%$$%poolMethod == "mice" && a%$$%toCall == "glm")  {
             retList <- ana
        }  else  {
             if( length(table(ana[,a%$$%imp])) > 1 ) {
                 retList <- jk2.pool ( datLong = ana, allNam=a[a%$$%allNam], forceSingularityTreatment = a%$$%forceSingularityTreatment, modus=mod)
             }  else  {
                 retList <- ana[,-match(c(a%$$%nest, a%$$%imp), colnames(ana))]
             }
        }
    ### p-Werte ergaenzen ... geschieht erst nach dem poolen der Standardfehler
        retList <- addSig (dat = retList, allNam = a[a%$$%allNam])
        return(retList)}

checkEngine <- function ( a) {
          if (a%$$%engine == "BIFIEsurvey") {
              if (!is.null(a%$$%nest) ) {
                  message("Engine 'BIFIEsurvey' currently does not work for nested imputation. Set 'engine' to 'survey'.")
                  a$engine <- "survey"
              }
              if ( length(unlist(lapply(c("glm", "quantile"), FUN = function ( y ) {unlist(lapply(c(a%$$%modus, a%$$%toCall), FUN = function (w) { grep(y, w)}))}))) > 0 ) {
                  message("Engine 'BIFIEsurvey' currently does not work for regression models and quantiles. Set 'engine' to 'survey'.")
                  a$engine <- "survey"
              }
              if ( !a%$$%type %in% c("JK2", "JK1", "NONE") ) {
                  message("Engine 'BIFIEsurvey' currently only works for jackknife 1 and jackknife 2. Set 'engine' to 'survey' due to type = '",a%$$%type,"'.")
                  a$engine <- "survey"
              }
              if ( !is.null(a%$$%adjust)) {
                   message("Engine 'BIFIEsurvey' currently does not work for adjusted means. Set 'engine' to 'survey'.")
                   a$engine <- "survey"
              }                                                                 ### bifie survey verbieten, wenn frequency table mit chi.square = TRUE und group.differences = TRUE
              if ( a%$$%toCall == "table" && !is.null(a%$$%group.differences.by) ) {
                   message("Engine 'BIFIEsurvey' currently does not work for chi.square tests of group differences of frequency tables. Set 'engine' to 'survey'.")
                   a$engine <- "survey"
              }
          }
          return(a%$$%engine)}

### Hilfsfunktion fuer eatRep()
checkWecForUV <- function(dat, allNam) {
          if(!is.null(attr(dat, "wrapperForWec"))) {
             auchUV <- allNam[["independent"]]
          }  else  {
             auchUV <- NULL
          }
          return(auchUV)}

checkJK.arguments <- function(a) {
          if ( a%$$%type %in% c("JK1", "JK2")) {
               if ( is.null(a%$$%repWgt)) {
                   if (is.null(a%$$%PSU) ) {
                       stop("For type = '",a%$$%type,"', 'PSU' must be defined if no replicate weights are specified (i.e., if 'repWgt' is NULL).")
                   }
                   if (a%$$%type == "JK2" && is.null(a%$$%repInd) ) {
                       stop("For type = '",a%$$%type,"', 'repInd' must be defined if no replicate weights are specified (i.e., if 'repWgt' is NULL).")
                   }
               }
          }}

checkGroupVars <- function ( datL, allNam, auchUV) {
    ### Personen mit Gewicht = 0 ausschliessen, damit sie nicht bei Ncases mitgezaehlt werden
    ### Aber Achtung: Gewicht = NA nicht mit ausschliessen, denn das ist kritisch und soll spaeter ggf. fehlermeldung werfen
          if(!is.null(allNam[["wgt"]])) {
             w0 <- which(datL[,allNam[["wgt"]]] == 0)
             if(length(w0) >0){
                nID  <- unique(datL[w0,allNam[["ID"]]])
                message("Remove ",length(nID), " students with zero weights to avoid counting them when determining sample size.")
                datL <- datL[-w0,]
             }
          }
          if(!is.null(allNam[["group"]]) || !is.null(auchUV) ) {
             chk <- lapply(allNam[["group"]], FUN = function ( v ) { if ( !inherits(datL[,v],  c("factor", "character", "logical", "integer"))) {stop(paste0("Grouping variable '",v,"' must be of class 'factor', 'character', 'logical', or 'integer'.\n"))} })
             for ( gg in c(allNam[["group"]], auchUV) ) {                       ### levels der Gruppen duerfen keine "." oder "_" enthalten, falls cross differences berechnet werden sollen
    ### personen muessen in gruppierungsvariable genestet sein, aber fuer jede Imputation und jedes Nest separat! (es sei denn, es gibt keine Imputationen oder die Gruppenvariable ist nicht imputiert)
                   if ( length(which(is.na(datL[,gg])))>0) { stop("Grouping variable '",gg,"' contains ",length(which(is.na(datL[,gg])))," missing values.")}
                   isImp<- TRUE                                                 ### initialisieren: wir gehen davon aus, dass die Variable imputiert ist
                   if ( is.null(allNam[["nest"]]) && is.null(allNam[["imp"]]) ) {
                        isImp <- FALSE                                          ### rauskriegen, ob fuer jede Imputation
                   }  else  {                                                   ### die genestetheit ueberprueft werden muss oder ob das fuer alle
                        imps  <- c(allNam[["nest"]], allNam[["imp"]])           ### zusammen geht
                        if ( length(imps) == 1) {
                             datL[,"g"] <- datL[,imps]
                        }  else  {
                             txt   <- paste0("datL |> tidyr::unite(\"g\", c(", paste(match(c(allNam[["nest"]],allNam[["imp"]]), colnames(datL)),collapse=", "),"), remove=FALSE)")
                             datL  <- eval(parse(text=txt))
                        }
                        isUniq<- eatGADS::checkUniqueness2(datL, varName=gg, idVar=allNam[["ID"]], impVar="g")
                        if(is.na(isUniq)) {
                            warning("'eatGADS::checkUniqueness2' returns NA. TRUE/FALSE is expected.", immediate. = TRUE)
                            isImp <- TRUE
                        }  else  {
                            if(isUniq) { isImp <- FALSE}
                        }
                   }
    ### wenn die Gruppierungsvariable nicht imputiert ist, genuegt es, den check nur fuer die erste Imputation auszufuehren
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
                   if (isFALSE(chk2)) { warning("Grouping variable '",gg,"' is not nested within persons (variable '",allNam[["ID"]],"').", immediate. = TRUE) }
    ### Umlaute-Bug aus dem LV 2012: workaround
                   if(inherits ( datL[,gg], c("factor", "character"))) {
                      if(inherits(try(datL[,gg] <- eatTools::cleanifyString(x = datL[,gg])  ),"try-error"))  {datL[,gg] <- eatTools::cleanifyString(x = datL[,gg], oldEncoding = "latin1")}
                   }
             }
          }
          if(!is.null(allNam[["group"]]) | !is.null(allNam[["independent"]]) ) {
             for ( gg in c(allNam[["group"]], allNam[["independent"]]) ) {
                 if (inherits ( datL[,gg], "factor")) {                         ### ausserdem duerfen fuer Gruppierungs- und unabhaengig Variablen keine factor levels ohne Beobachtungen drinsein
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

checkForAdjustmentAndLmer <- function(datL, a, groupWasNULL) {
    ### adjustierungsvariablen und lmer-Praediktoren fuer multilevel analysen in bifie survey duerfen nur numerisch oder dichotom sein und sollten keine missings enthalten
          if(!is.null(a%$$%adjust) || !is.null(a%$$%formula.random) || !is.null(a%$$%formula.fixed) ) {
             vars<- rbind(data.frame (type = rep("adjust", length(a%$$%adjust)),vars = a%$$%adjust , stringsAsFactors = FALSE),
                          data.frame (type = rep("fixed", length(all.vars(a%$$%formula.fixed))),vars = all.vars(a%$$%formula.fixed) , stringsAsFactors = FALSE),
                          data.frame (type = rep("random", length(all.vars(a%$$%formula.random))),vars = all.vars(a%$$%formula.random) , stringsAsFactors = FALSE))
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
                       if (vars[g,"type"] == "adjust") {                        ### untere Zeile: Levels der Adjustierungsvariablen duerfen keine Leerzeichen, Doppelpunkte, Klammer auf oder Klammer zu haben, da model.matrix die als colnames verwendet und dort sind sie in data.frames nicht erlaubt
                            datL[,vars[g,"vars"]] <- gsub("\\(|\\)|:| |-|,", ".", as.character(datL[,vars[g,"vars"]]))
                            newFr <- model.matrix( as.formula (paste("~",vars[g,"vars"],sep="")), data = datL)[,-1,drop=FALSE]
                            message(paste("    Adjusting variable '",vars[g,"vars"],"' of class '",class(datL[,vars[g,"vars"]]),"' was converted to ",ncol(newFr)," indicator(s) with name(s) '",paste(colnames(newFr), collapse= "', '"), "'.",sep=""))
                            datL  <- data.frame(datL, newFr, stringsAsFactors=FALSE)
                            a$adjust <- c(setdiff(a%$$%adjust,vars[g,"vars"]), colnames(newFr))
                       }
                       if (vars[g,"type"] == "fixed") {
                            if ( !vars[g,"vars"] %in% extractFactorVarsFromFormula(a%$$%formula.fixed) && vars[g,"vars"] %in% all.vars(a%$$%formula.fixed) ) {
                                 stop(paste0("Variable '",vars[g,"vars"],"' of class '",class(datL[,vars[g,"vars"]]), "' should be specified with 'as.factor(",vars[g,"vars"],")' in the formula.fixed argument.\n  If '",vars[g,"vars"],"' should be modeled as numeric, please change variable class of '",vars[g,"vars"],"' in the data into numeric."))
                            }
                       }
                       if (vars[g,"type"] == "random") {
                            if ( !vars[g,"vars"] %in% extractFactorVarsFromFormula(a%$$%formula.random) && vars[g,"vars"] %in% all.vars(a%$$%formula.random) ) {
                                 warning(paste0("Variable '",vars[g,"vars"],"' of class '",class(datL[,vars[g,"vars"]]), "' should be specified with 'as.factor(",vars[g,"vars"],")' in the formula.random argument.\n  If '",vars[g,"vars"],"' should be modeled as numeric, please change variable class of '",vars[g,"vars"],"' in the data into numeric."), immediate. = TRUE)
                            }
                       }
                 }
             }
    ### wenn es adjustierungsvariablen gibt, darf groups nicht NULL gewesen sein, also nicht "wholeGroup" sein ... es darf aber (zumindest zunaechst einmal) nur eine Gruppierungsvariable geben
             if ( groupWasNULL && !is.null(a%$$%adjust)) {stop("When adjusted variables are defined, argument 'groups' must not be NULL.")}
             if ( length(a%$$%group)>1  && !is.null(a%$$%adjust) ) {stop("When adjusted variables are defined, to date, only one grouping variable is allowed.")}
    ### wenn es adjustierungsvariablen gibt, koennen cross-level differences vorerst nur mit methode 'old' bestimmt werden
    ### dieser check findet in der Funktion repMean statt ... aber: wenn es adjustierungsvariablen gibt, kann nicht BIFIEsurvey benutzt werden
          }
          a[["datL"]] <- datL
          return(a)}
          
checkNameConvention <- function( allNam) {
          na    <- c("N_weightedValid", "N_weighted",  "wgtOne", "wgtOne2","le", "variable", "est", "se")
          naGr  <- c("wholePop", "group", "depVar", "modus", "parameter", "coefficient", "value", "linkErr", "comparison", "sum", "trendvariable", "g", "le", "splitVar", "rowNr", "variable", "Freq", "id", "unit_1", "unit_2", "comb.group", "row", "coef", "label1", "label2", "hierarchy.level")
          naInd <- c("(Intercept)", "Ncases", "Nvalid", "R2",  "R2nagel", "linkErr", "variable")
          naGr1 <- which ( allNam[["group"]] %in% naGr )                        ### hier kuenftig besser: "verbotene" Variablennamen sollen automatisch umbenannt werden!
          if(length(naGr1)>0)  {stop(paste0("Following name(s) of grouping variables in data set are forbidden due to danger of confusion with result structure:\n     '", paste(allNam[["group"]][naGr1], collapse="', '"), "'\n  Please rename these variable(s) in the data set.\n"))}
          naInd1<- which ( allNam[["independent"]] %in% naInd )
          if(length(naInd1)>0)  {stop(paste0("Following name(s) of independent variables in data set are forbidden due to danger of confusion with result structure:\n     '", paste(allNam[["independent"]][naInd1], collapse="', '"), "'\n  Please rename these variable(s) in the data set.\n"))}
          na2   <- which ( unlist(allNam) %in% na )
          if(length(na2)>0)  {stop(paste0("Following variable name(s) in data set are forbidden due to danger of confusion with result structure:\n     '", paste(unlist(allNam)[na2], collapse="', '"), "'\n  Please rename these variable(s) in the data set.\n"))} }

setCrossDifferences <- function (a) {
          if ( !is.list(a%$$%cross.differences)) {
               if ( !is.logical (a%$$%cross.differences)) {
                    message("Argument 'cross.differences' must be either logical or a list. Set 'cross.differences' to FALSE.")
                    a$cross.differences <- FALSE
               }
               if ( isTRUE(a%$$%cross.differences)) {
                    if ( "wholeGroup" %in% (a%$$%group) ) {
                          message("No groups defined. Set 'cross.differences' to FALSE.")
                          a$cross.differences <- FALSE
                    }  else  {
                          if ( length(a%$$%group.splits)==1) {
                                message("Argument 'group.splits' was set to 1. No cross-level differences can be computed. Set 'cross.differences' to FALSE.")
                                a$cross.differences <- FALSE
                          }  else  {
                                a$cross.differences <- combinat::combn(x=a%$$%group.splits,m=2, simplify = FALSE)
                          }
                    }
               }
          }  else  {                                                            ### checks muessen nur gemacht werden, wenn die Liste nicht automatisch erzeugt wurde
               if(!all(unlist(lapply(a%$$%cross.differences, length)) == 2)) {stop("Each element in 'cross.differences' must be a vector of length 2.\n")}
               if(!all(unlist(lapply(a%$$%cross.differences, inherits, what = c("integer","numeric"))))) {stop("Each element in 'cross.differences' must be a numerical vector.\n")}
               if(!all(unlist(lapply(a%$$%cross.differences, FUN = function ( x ) { length(unique(x))})) ==2 ) ) {stop("Each element in 'cross.differences' must be a vector of 2 different numbers.\n")}
               vals <- unique(unlist(a%$$%cross.differences))
               if(!all(vals %in% (a%$$%group.splits))) {stop("All numerical values in 'cross.differences' must be included in 'group.splits'.\n")}
               a$cross.differences <- lapply( a%$$%cross.differences, sort )
               if ( length(a%$$%cross.differences) != length(unique(a%$$%cross.differences)) ) {
                    message("Some comparisons in 'cross.differences' are identical. Duplicated comparisons will be removed.")
                    a$cross.differences <- unique(a%$$%cross.differences)
               }
          }
          a[["allNam"]] <- c(a[["allNam"]], "cross.differences")
          return(a)}

### Achtung: wenn keine Gruppen und/oder Nests und/oder Imputationen spezifiziert sind, erzeuge Variablen mit Werten gleich 1, damit by() funktioniert!
createLoopStructure <- function (a) {
          if( is.null(a%$$%imp) )  { a$datL[,"imp"] <- 1; a$imp <- "imp" } else { stopifnot(length(a%$$%imp) == 1 ); a$datL[,a$imp] <- as.character(a$datL[,a$imp])}
          if( is.null(a%$$%wgt) )  { a$datL[,"wgtOne"] <- 1; a$wgt <- "wgtOne" } else {
              stopifnot(length(a%$$%wgt) == 1 )
              if ( !inherits(a[["datL"]][,a%$$%wgt], c("numeric", "integer")) ) { stop ( paste("Error: 'wgt' variable '",a%$$%wgt,"' of class '",paste(class(a[["datL"]][,a%$$%wgt]),collapse = "', '"),"' has to be numeric.\n",sep="")) }
              isMis <- which(is.na(a[["datL"]][,a%$$%wgt]))
              isZero<- which ( a[["datL"]][,a%$$%wgt] == 0 )
              if(length(isMis)>0) { stop (paste ( "Error: Found ",length(isMis)," missing values in the weight variable '",a%$$%wgt,"'.\n",sep="")) }
              if(length(isZero)>0) { warning( "Found ",length(isZero)," zero weights in the weight variable '",a%$$%wgt,"'.", immediate. = TRUE) }
          }
          if( is.null(a%$$%L2wgt) )  { a[["datL"]][,"wgtOne2"] <- 1; a$L2wgt <- "wgtOne2" }
          if(!is.null(a%$$%nest))  {
              stopifnot(length(a%$$%nest) == 1 )
              a[["datL"]][,a%$$%nest] <- as.character(a[["datL"]][,a%$$%nest])
              if(a%$$%verbose){cat(paste("\nAssume nested structure with ", length(table(a[["datL"]][,a%$$%nest]))," nests and ",length(table(a[["datL"]][,a%$$%imp]))," imputations in each nest. This will result in ",length(table(a[["datL"]][,a%$$%nest]))," x ",length(table(a[["datL"]][,a%$$%imp]))," = ",length(table(a[["datL"]][,a%$$%nest]))*length(table(a[["datL"]][,a%$$%imp]))," imputation replicates.\n",sep=""))}
          }  else  { if(a%$$%verbose){cat("\nAssume unnested structure with ",length(table(a[["datL"]][,a%$$%imp]))," imputations.\n",sep="")}}
          if( is.null(a%$$%nest) ) { a$datL[,"nest"]  <- 1; a$nest  <- "nest" }
          if(!is.null(a%$$%group)) {                                            ### untere Zeile: das, damit leere Gruppen nicht ueber by() mit geschleift werden, wie es passiert, wenn Gruppen als Faktoren definiert sind
              for ( jj in a%$$%group )  { a$datL[,jj] <- as.character(a$datL[,jj]) }
          }
          return(a)}

### replicates zuweisen bzw. erzeugen
assignReplicates <- function ( a) {
    ### wenn Replicates bereits uebergeben, muss PSU und repInd NULL sein
          if(!is.null(a%$$%repWgt) ) {
              if ( !is.null(a%$$%PSU) | !is.null(a%$$%repInd) ) {
                    warning("Arguments 'PSU' and 'repInd' are expected to be NULL if replicate weights are already defined (via 'repWgt').\n    'PSU' and 'repInd' will be ignored.", immediate. = TRUE)
              }
          }
    ### replicates erzeugen (nur einmal fuer alle Analysen, und nur wenn 'repWgt' NULL ist)
          if(!is.null(a%$$%repWgt))  {
              repA <- data.frame ( a[["datL"]][,a%$$%ID, drop=FALSE], a[["datL"]][,a%$$%repWgt ])
              repA <- repA[!duplicated(repA[,a%$$%ID]),]
          }  else  {
              if(!is.null(a%$$%PSU) && a%$$%engine=="survey" )  {
                  repW <- a[["datL"]][!duplicated(a[["datL"]][,a%$$%ID]),]
                  repA <- generate.replicates(dat = repW, a=a)
              }  else  { repA <- NULL}
          }
          return(repA)}

generate.replicates <- function ( dat, a)   {
          for ( i in names(a)) { assign(i, a[[i]]) }
          if(type %in% c("JK2", "BRR")) { stopifnot(length(PSU) == 1 & length(repInd) == 1 ) }
          if(type  == "JK1" ) { if(!is.null(repInd))  {
             cat("'repInd' is ignored for 'type = JK1'.\n")
             a["repInd"] <- list(NULL)
          }  }
          vars        <- unlist(a[allNam][c("ID", "wgt", "PSU", "repInd")])
          dat.i       <- dat[,vars]
          if(type %in% c("JK2", "BRR")) { if( !all( names(table(dat.i[,repInd])) == c(0,1)) ) {stop("Only 0 and 1 are allowed for 'repInd' variable.\n")} }
          zonen       <- names(table(as.character(dat.i[,PSU]) ) )
          if ( verbose) { cat(paste("Create ",length(zonen)," replicate weights according to ",type," procedure.\n",sep=""))}
          if ( progress && nrow(dat)>5000 & length(zonen) > 50 ) {              ### progress bar wird nur angezeigt, wenn es auch nennenswerte Zeit dauert
               pb     <- progress::progress_bar$new( format = "         replicates [:bar] :percent in :elapsed", incomplete = " ", total = length(zonen), clear = FALSE, width= 75, show_after = 0.01)
          }
          missings    <- sapply(dat.i, FUN = function (ii) {length(which(is.na(ii)))})
          if(!all(missings == 0)) {
              mis.vars <- paste(names(missings)[which(missings != 0)], collapse = ", ")
              stop(paste("Found missing value(s) in variable(s) ", mis.vars,".\n",sep=""))
          }                                                                     ### untere Zeile: in JK-2 werden die gewichte nur in der entsprechenden PSU FUER EIN REPLICATE geaendert
          reps <- data.frame ( lapply(zonen , FUN = function(ii) {              ### in BRR werden die Gewichte in allen Zonen FUER EIN REPLICATE geaendert!
                  if ( progress && nrow(dat)>5000 & length(zonen) > 50 ) { pb$tick() }
                  rep.ii <- dat.i[,wgt]                                         ### in JK-1 werden die Gewichte nur in der entsprechenden PSU FUER ALLE REPLICATES geaendert
                  if(type == "JK2")  { rep.ii[dat.i[,PSU] == ii ] <- ifelse(dat.i[ dat.i[,PSU] == ii ,repInd] == 1, 0, 2 * rep.ii[dat.i[,PSU] == ii ] ) }
                  if(type == "BRR")  { rep.ii <- ifelse(dat.i[ ,repInd] == 1, 0, 2 * rep.ii ) }
                  if(type == "JK1")  {
                     rep.ii[ which ( dat.i[,PSU] == ii) ] <- 0
                     rep.ii[ which ( dat.i[,PSU] != ii) ] <- rep.ii[ which ( dat.i[,PSU] != ii) ] *  ( sum(dat.i[,wgt]) / sum (rep.ii))
                  }
                  return(rep.ii) }), stringsAsFactors = FALSE)
          colnames(reps) <- paste(wgt, 1:ncol(reps), sep="_")
          ret            <- data.frame(dat.i[,ID,drop=FALSE], reps, stringsAsFactors = FALSE)
          attr(ret, "n.replicates") <- length(zonen)
          return(ret) }

checkImpNest <- function (toAppl, gr, a) {
          for ( i in names(a)) { assign(i, a[[i]]) }
          if(isTRUE(doCheck)) {                                                 ### wenn check=TRUE, werden fuer jede Analyse des splitters
             if ( length( toAppl[[gr]] ) > 1) {                                 ### eine Reihe von checks durchgefuehrt
                  crsTab <- table(datL[,toAppl[[gr]]])
                  if ( length(which(crsTab < 10 )) > 0 ) {
                       warning("Small number of observations in some combinations of grouping variables:\n   Recommend to remove these group(s).\n", eatTools::print_and_capture(crsTab, 3) , immediate. = TRUE)
                  }
             }
             impNes<- by(data = datL, INDICES = datL[, c(nest, toAppl[[gr]]) ], FUN = function ( x ) { length(table(as.character(x[,imp])))}, simplify = FALSE)
             laenge<- which(sapply(impNes, length) == 0)
             if ( length(laenge ) > 0 ) {
                  warning(paste0(length(laenge), " combination(s) of groups without any observations. Analysis might crash."), immediate. = TRUE)
             }
    ### check: gleich viele nests/imputationen je kombination von gruppierungsvariablen?
             chk1  <- nestsImpsPerGroupComb(datL=datL, allNam=a[allNam], toAppl=toAppl, gr=gr)
    ### check: gleichviele PSUs je nest?
             if(!is.null(PSU))  {
                  psuNes<- table ( by(data = datL, INDICES = datL[,nest], FUN = function ( x ) { length(table(as.character(x[,PSU])))}, simplify = FALSE) )
                  if(length(psuNes) != 1 ) {warning("Number of PSUs differ across nests!\n", eatTools::print_and_capture(psuNes, 3), immediate. = TRUE)}
             }
    ### check: sind fuer jede Gruppe alle Faktorstufen in allen nests und allen imputationen vorhanden? z.B. nicht in einer Imputation nur Jungen
             impNes<- by(data = datL, INDICES = datL[, c(nest, imp) ], FUN = checkNests, allNam=a[allNam], toAppl=toAppl, gr=gr, simplify = FALSE)
             impNes<- data.frame ( do.call("rbind", lapply(impNes, FUN = function ( x ) { unlist(lapply(x[["ret"]], FUN = length)) })) )
             if ( !all ( sapply(impNes, FUN = function ( x ) { length(table(x)) } ) == 1) ) { warning("Number of units in at least one group differs across imputations!", immediate. = TRUE)}
    ### Achtung!! jetzt der check, der in der alten Version ueber 'checkData' gemacht wurde!
             chk2  <- do.call("rbind", by(data = datL, INDICES = datL[,c( group, nest, imp)], FUN = checkData, a=a))
          }}

### Hilfsfunktion fuer checkImpNest
nestsImpsPerGroupComb <- function(datL, allNam, toAppl, gr) {
         impNes<- table(datL[, c(allNam[["nest"]], allNam[["imp"]], toAppl[[gr]]) ])
         inDF  <- as.data.frame(impNes)
         if ( length(toAppl[[gr]]) >0) {ind <- inDF[,toAppl[[gr]]]} else {ind <- rep(1, nrow(inDF))}
         grpVar<- by(data=inDF, INDICES = ind, FUN = function (grp) {
                  ch1 <- by(data = grp, INDICES = grp[,allNam[["nest"]]], FUN = function (i) {unique(i[,allNam[["imp"]]])})
                  if(length(ch1)>1) {
                      cmb <- combinat::combn(x=1:length(ch1), m=2, simplify=FALSE)# gleichviele imputationen je nest?
                      ch1 <- all(unlist(lapply(cmb, FUN = function (j) {isTRUE(all.equal(ch1[[j[1]]], ch1[[j[2]]]))})))
                  }  else {
                      ch1 <- TRUE
                  }
                  ch2 <- all(grp[,"Freq"]>0)
                  ch3 <- length(unique(grp[,"Freq"])) == 1
                  return(c(ch1, ch2, ch3))})
         if ( any(impNes==0) || any( unlist(grpVar) == FALSE) ) {warning("Number of imputations differ across nests and/or groups!\n", eatTools::print_and_capture(impNes, 3), immediate. = TRUE)}}

prepExpecVal <- function (a) {
          if(a%$$%toCall=="table") {
             misInd <- which(is.na(a$datL[,a%$$%dependent]))
             if(isTRUE(a%$$%separate.missing.indicator)) {
                if(length(misInd)>0) { a$datL[misInd,a%$$%dependent] <- "<NA>"}
             }  else {
                if(length(misInd)>0) {
                   warning("No seperate missing categorie was chosen. ", length(misInd), " missings were found anyhow for ",a%$$%dependent,". Missings will be deleted from the data.", immediate. = TRUE)
                   if(length(misInd) == nrow(a%$$%datL)) {stop()}
                   a$datL <- a$datL[-misInd,]
                }
             }
             a$expected.values <- sort(unique(c(a%$$%expected.values, names(table(a$datL[,a%$$%dependent])))))
          }
          return(a)}

createInfoString <- function ( ai, toCall, gr, toAppl) {
          nr <- ifelse(is.null(ai), 1, match(gr, names(toAppl)))
          st <- paste0(paste(rep(" ", times = 8-nchar(toCall)),collapse=""), toCall, " analysis ",nr)
          return(st)}
          
keepNonNULL <- function(list1, list2){
          comm <- intersect(names(list1), names(list2))
          vgl  <- data.frame ( name = comm, list1 = unlist(lapply(comm, FUN = function (co) {is.null(list1[[co]])})), list2 = unlist(lapply(comm, FUN = function (co) {is.null(list2[[co]])})), stringsAsFactors = FALSE)
          weg1 <- intersect(which(vgl[,"list1"]), which(!vgl[,"list2"]))
          if(length(weg1)>0) {list1 <- list1[-eatTools::whereAre(vgl[weg1,"name"], names(list1), verbose=FALSE)]}
          weg2 <- intersect(which(vgl[,"list2"]), which(!vgl[,"list1"]))
          if(length(weg2)>0) {list2 <- list2[-eatTools::whereAre(vgl[weg2,"name"], names(list2), verbose=FALSE)]}
          if(length(intersect(names(list1), names(list2))) > 0 ) {
              ret  <- c(list1[-eatTools::whereAre(names(list2), names(list1), verbose=FALSE)], list2)
          }  else  {
              ret  <- c(list1, list2)
          }
          return(ret)}



computeTrend <- function(jk2, tv, repFunOut, fun) {
        jk2_all <- do.call("rbind", jk2)                                        ### bind yearwise results
        jk2_bind<- jk2_all
   ### checks: die selben Zeilen in allen Jahren?? (fuer GLMs insbesondere testen!)
        jk2_bind<- check1(jk2=jk2, jk2_bind=jk2_bind, jk2_all=jk2_all, tv=tv)
   ### special for mean: select only mean and sd, reshape le
        if(identical(fun, "mean")) {
            jk2_bind<- jk2_bind[jk2_bind[["parameter"]] %in% c("mean", "sd"), ]
        }
        if ( fun == "glm" ) {
            jk2_bind<- jk2_bind[which(!jk2_bind[,"parameter"] %in% c("Nvalid", "Ncases", "R2", "R2nagel")),]
        }
        jk2_bind<- jk2_bind[jk2_bind[["coefficient"]] %in% c("est", "se"), ]    ### drop significance
        jk2_wide<- reshape2::dcast ( jk2_bind, as.formula ( paste ( " ... ~ coefficient + ", tv ,sep="") ) )
        lev     <- unique(jk2_bind[[tv]])                                       ### calculate trend
        le      <- check2(repFunOut=repFunOut, jk2=jk2_bind, fun=fun)
   ### calculate all trends!
        vgl <- combinat::combn(names(jk2),2, simplify=FALSE)                    ### untere Zeile: Zusaetze fuer alle trendvergleiche (1 vs. 2; 2 vs. 3; 1 vs. 3) machen
        adds<- lapply( 1:length(vgl), FUN = function ( comp) {                  ### untere Zeile: merge linking errors ... muss angepasst werden!
               jk2_wide[,paste0("est_trend_",vgl[[comp]][1],".vs.",vgl[[comp]][2])] <- jk2_wide[,paste("est_",vgl[[comp]][2],sep="")] - jk2_wide[,paste("est_",vgl[[comp]][1],sep="")]
   ### merge linking errors ... aus 'le' die relevanten Jahre auswaehlen: wenn der user die in umgekehrter Reihenfolge, also 2015 in die erste, und 2010 in die zweite Spalte geschrieben hat, muss das jetzt homogenisiert werden
               years.c <- unlist(apply(le, MARGIN = 1, FUN = function (zeile) {paste(c(zeile[["trendLevel1"]], zeile[["trendLevel2"]]), collapse="_")}))
               le.years<- le[which(years.c == paste(sort(vgl[[comp]]), collapse="_")),]
               jk2_wide<- merge(jk2_wide, le.years[,-na.omit(match(c("trendLevel1", "trendLevel2", "domain"), colnames(le.years)))], by = c("parameter", "depVar"), all = TRUE)
   ### check for missing LEs
               miss    <- which(is.na(jk2_wide[,"le"]))
               if ( length(miss)>0){
                    warning(paste0("Found ",length(miss)," missing linking errors for dependent variable '",unique(jk2_wide[,"depVar"]),"' and parameter(s) '",paste(unique(jk2_wide[which(is.na(jk2_wide[,"le"])),"parameter"])),"'. Assume linking error of 0 for these cases."))
                    jk2_wide[which(is.na(jk2_wide[,"le"])),"le"] <- 0
               }
   ### calculate trend SEs and p values
               jk2_wide[,paste0("se_trend_", vgl[[comp]][1],".vs.",vgl[[comp]][2])] <- sqrt(jk2_wide[, paste("se_",vgl[[comp]][1],sep="")]^2 + jk2_wide[, paste("se_",vgl[[comp]][2],sep="")]^2 + jk2_wide[, "le"]^2)
               jk2_wide[,paste0("sig_trend_", vgl[[comp]][1],".vs.",vgl[[comp]][2])]<- 2*pnorm(abs(jk2_wide[, paste0("est_trend_",vgl[[comp]][1],".vs.",vgl[[comp]][2])]/jk2_wide[, paste0("se_trend_", vgl[[comp]][1],".vs.",vgl[[comp]][2])]), lower.tail=FALSE)
   ### Effect size for means (only for jk2.mean)
               existSD <- "sd" %in% jk2_bind[,"parameter"]                      ### nur wenn standardabweichungen drin stehen, koennen effektstaerken berechnet werden ... Message wird nur fuer ersten Schleifendurchlauf ausgegeben (redundanz vermeiden)
               es  <- character(0)
               if(fun == "mean" && !existSD && comp == 1) {message("Cannot find standard deviations in output. Skip computation of effect sizes.")}
               if (  fun == "mean" && existSD) {                                ### not for groupDiffs as no SD is provided by eatRep, split up data frame and rbind later
                   jk2_wide[, paste0("es_trend_", vgl[[comp]][1],".vs.",vgl[[comp]][2])]  <- NA
                   jk2_wideS<- jk2_wide[!jk2_wide[, "comparison"] %in% c("crossDiff_of_groupDiff", "groupDiff") | is.na(jk2_wide[, "comparison"]), ]
                   tabs     <- table(jk2_wideS[,c("parameter", "group")])
                   if ( !all ( tabs == 1) ) {                                   ### checks SW: jeder Eintrag in "group" darf nur zweimal vorkommen ... sonst Effektstaerkebestimmug ueberspringen
                        message("Cannot find standard deviations of following groups: \n",print_and_capture(tabs, 5),"\nSkip computation of effect sizes.")
                   }  else  {                                                   ### sortieren nach "group" und "parameter"; siehe Mail an Benjamin, 11.06.2019
                        jk2_wideS<- jk2_wideS[with(jk2_wideS, order(parameter, group)),]
                        pooledSD <- sqrt(0.5 * (jk2_wideS[jk2_wideS[, "parameter"] == "sd", paste("est_",vgl[[comp]][1], sep="")]^2 + jk2_wideS[jk2_wideS[, "parameter"] == "sd", paste("est_",vgl[[comp]][2], sep="")]^2))
                        jk2_wideS[jk2_wideS[, "parameter"] == "mean", paste0("es_trend_", vgl[[comp]][1],".vs.",vgl[[comp]][2])]  <- jk2_wideS[jk2_wideS[, "parameter"] == "mean", paste0("est_trend_", vgl[[comp]][1],".vs.",vgl[[comp]][2])] / pooledSD
                        es       <- c(es, paste0("es_trend_", vgl[[comp]][1],".vs.",vgl[[comp]][2]))
                        jk2_wide <- unique(rbind(jk2_wide, jk2_wideS))          ### add es to reshaping!
                   }
               }
   ### reshape back to standard format
               jk2_add <- reshape2::melt ( jk2_wide, measure.vars = c(paste("est", unique(unlist(vgl)), sep="_"),  paste("se", unique(unlist(vgl)), sep="_"),  paste0("est_trend_",vgl[[comp]][1],".vs.",vgl[[comp]][2]), paste0("se_trend_",vgl[[comp]][1],".vs.",vgl[[comp]][2]), paste0("sig_trend_",vgl[[comp]][1],".vs.",vgl[[comp]][2]), es), na.rm = FALSE)
               jk2_add <- jk2_add[!is.na(jk2_add[, "value"]), ]                  ### drop missing rows (e.g. for SD, as no es can be calc)
   ### split up variable column into coef and trend variable ... an ERSTEM Unterstrich splitten
               zusatz  <- data.frame ( eatTools::halveString (string = as.character(jk2_add[,"variable"]), pattern = "_", first = TRUE ), stringsAsFactors = FALSE)
               colnames(zusatz) <- c("coefficient", tv)
               jk2_add <- data.frame ( jk2_add[,-match("variable", colnames(jk2_add))], zusatz, stringsAsFactors = FALSE)
   ### return additional results (drops le)
               return(jk2_add[, names(jk2_all)])})
   ### Originaloutput und Zusatz untereinanderbinden, nur das unique behalten
        jk2 <- unique(rbind(jk2_all, do.call("rbind", adds)))
        return(jk2) }

check1 <- function(jk2, jk2_bind, jk2_all, tv) {
    ### alles gegen alles vergleichen
       vgl <- combinat::combn(1:length(jk2),2, simplify=FALSE)
       chks<- sapply(vgl, FUN = function ( x ) { length(unique(jk2[[x[1]]]$group)) != length(unique(jk2[[x[2]]]$group)) || suppressWarnings(!all(sort(unique(jk2[[x[1]]]$group)) == sort(unique(jk2[[x[1]]]$group))))})
       if(sum(chks)>0) {                                                        ### wenn wenigstens ein check TRUE ist, soll das gemacht werden
           weg <- unique(unlist(lapply(vgl, FUN = function ( v ) {              ### remove rows which are not in all years (trend only for rows which are in both data sets)
                  uneven1 <- setdiff(unique(jk2[[v[1]]]$group), unique(jk2[[v[2]]]$group))
                  uneven2 <- setdiff(unique(jk2[[v[2]]]$group), unique(jk2[[v[1]]]$group))
                  if(length(uneven1) > 0) {cat("Categories", uneven1, "are missing in year", jk2[[2]][1, tv] ,". \n")}
                  if(length(uneven2) > 0) {cat("Categories", uneven2, "are missing in year", jk2[[1]][1, tv] ,". \n")}
                  weg1    <- which( jk2_all$group %in% c(uneven1, uneven2))
                  return(weg1)})))
           if ( length(weg)>0) {
                jk2_bind <- jk2_all[-weg,]
           }
       }
       return(jk2_bind)}

check2 <- function(repFunOut, jk2, fun){
    ### wenn linking error objekt urspruenglich kein data.frame, dann hier einfach durchschleifen
       wdf <- attr(repFunOut[["le"]], "linkingErrorFrame")                      ### was data.frame?
       if(is.null(wdf)) {
           return(repFunOut[["le"]])
       }  else  {
    ### notwendige Spalten
           allV  <- list(trendLevel1 = "trendLevel1", trendLevel2 = "trendLevel2", parameter = "parameter", linkingError="linkingError", depVar = "depVar")
           allN  <- lapply(allV, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = repFunOut[["le"]], variable=ii, warnIfMissing = TRUE)})
    ### AV in linking error objekt enthalten? Achtung! wenn repTable ueber repMean gewrappt wurde, heisst die AV anders. Dann (aber nur dann!) muss das le-Objekt angepasst werden
           le    <- repFunOut[["le"]]
           dv    <- unique(jk2[,"depVar"])
           stopifnot(length(dv)==1)
           if(!dv %in% le[,"depVar"]) {stop(paste0("Cannot found dependent variable '",dv,"' in 'depVar' column of linking error data.frame 'linkErr'."))}
#               if (fc != "repTable") {
#                   stop(paste0("Cannot found dependent variable '",allNam[["dependent"]],"' in 'depVar' column of linking error data.frame 'linkErr'."))
#               }  else  {
#                   lvls <- unique(le[,"depVar"])                             ### assumed level
#                   mat  <- lapply(lvls, FUN = function (l) {grep(l,allNam[["dependent"]])})
#                   mat1 <- which(mat == 1)
#                   if ( length(mat1) != 1) {stop(paste0("Cannot match dependent variable in data ('",allNam[["dependent"]],"') to the dependent variables in the 'depVar' column of linking error data.frame 'linkErr': '",paste(unique(le[,"depVar"]), collapse = "', '"),"'."))}
#                   le   <- le[which(le[,"depVar"] == lvls[mat1]),]
#                   prm  <- eatTools::removePattern( allNam[["dependent"]], pattern = lvls[mat1])
#                   if ( !prm %in% le[,"parameter"]) {stop(paste0("Cannot find parameter '",prm,"' in the 'parameter' column of linking error data.frame 'linkErr'"))}
#                   le   <- le[which(le[,"parameter"] == prm),]
#                   le[,"depVar"] <- allNam[["dependent"]]
#               }
#           }  else  {
#               le  <- le[which(le[,"depVar"] == allNam[["dependent"]]),]
#           }
           le  <- le[which(le[,"depVar"] == dv),]
    ### nur ein Kompetenzbereich?
           le1   <- le[,c("parameter", "trendLevel1", "trendLevel2")]
           if(nrow(le1) != nrow(unique(le1))) {stop("Linking error data.frame 'linkErr' is not unique in 'parameter', 'trendLevel1', 'trendLevel2'. Probable cause: 'linkErr' contains linking errors of multiple competence domains?")}
    ### gibt es fuer alle Parameter auch linkingfehler? Hier reicht eine warnung aus, da die Linkingfehler in der computeTrend-Funktion rangemergt werden. Da kann man sie dann auf 0 setzen
           add <- setdiff(unique(jk2[,"parameter"]), unique(le[,"parameter"]))
           if ( length(add)>0) { warning(paste0("No linking errors for parameters '",paste(add, collapse="', '"),"'. Linking errors for these parameters will be defaulted to 0."))}
    ### alle Kombinationen aus jahreszahlen vorhanden?
           years <- sort(unique (jk2[,repFunOut[["allNam"]][["trend"]]]))
           combs <- combinat::combn(x=years, m=2, simplify=FALSE)
           exist <- unique(plyr::alply(le, .margins = 1, .fun = function (zeile) {sort(c(zeile[["trendLevel1"]], zeile[["trendLevel2"]]))}))
           drin  <- combs %in% exist
           if(!all(drin)) {stop(paste0("Data contains combination of years '",paste(combs[!drin], collapse="', '"), "' which are not included in linking error data.frame 'linkErr'."))}
    ### linkingError spalte in le umbenennen
           colnames(le) <- car::recode(colnames(le), "'linkingError'='le'")
           return(le)
       }}


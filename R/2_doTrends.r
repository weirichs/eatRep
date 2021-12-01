computeTrend <- function(jk2, tv, le, fun) {
        jk2_all <- do.call("rbind", jk2)                                        ### bind yearwise results
        jk2_bind<- jk2_all
   ### checks: die selben Zeilen in allen Jahren?? (fuer GLMs insbesondere testen!)
        jk2_bind<- check1(jk2=jk2, jk2_bind=jk2_bind, jk2_all=jk2_all)
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
   ### calculate all trends!
        vgl <- combinat::combn(names(jk2),2, simplify=FALSE)                    ### untere Zeile: Zusaetze fuer alle trendvergleiche (1 vs. 2; 2 vs. 3; 1 vs. 3) machen
        adds<- lapply( 1:length(vgl), FUN = function ( comp) {                  ### untere Zeile: merge linking errors ... muss angepasst werden!
               jk2_wide[,paste0("est_trend_",vgl[[comp]][1],".vs.",vgl[[comp]][2])] <- jk2_wide[,paste("est_",vgl[[comp]][2],sep="")] - jk2_wide[,paste("est_",vgl[[comp]][1],sep="")]
               jk2_wide<- merge(jk2_wide, le, by = c("parameter", "depVar"), all = TRUE)
   ### calculate trend SEs and p values
               jk2_wide[,paste0("se_trend_", vgl[[comp]][1],".vs.",vgl[[comp]][2])] <- sqrt(jk2_wide[, paste("se_",vgl[[comp]][2],sep="")]^2 + jk2_wide[, paste("se_",vgl[[comp]][2],sep="")]^2 + jk2_wide[, "le"]^2)
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
               zusatz  <- data.frame ( halveString (string = as.character(jk2_add[,"variable"]), pattern = "_", first = TRUE ), stringsAsFactors = FALSE)
               colnames(zusatz) <- c("coefficient", tv)
               jk2_add <- data.frame ( jk2_add[,-match("variable", colnames(jk2_add))], zusatz, stringsAsFactors = FALSE)
   ### return additional results (drops le)
               return(jk2_add[, names(jk2_all)])})
   ### Originaloutput und Zusatz untereinanderbinden, nur das unique behalten
        jk2 <- unique(rbind(jk2_all, do.call("rbind", adds)))
        return(jk2) }
        
check1 <- function(jk2, jk2_bind, jk2_all) {
    ### alles gegen alles vergleichen
       vgl <- combinat::combn(1:length(jk2),2, simplify=FALSE)
       chks<- sapply(vgl, FUN = function ( x ) { length(unique(jk2[[x[1]]]$group)) != length(unique(jk2[[x[2]]]$group)) || suppressWarnings(!all(sort(unique(jk2[[x[1]]]$group)) == sort(unique(jk2[[x[1]]]$group))))})
       if(sum(chks)>0) {                                                        ### wenn wenigstens ein check TRUE ist, soll das gemacht werden
           weg <- unique(unlist(lapply(vgl, FUN = function ( v ) {                            ### remove rows which are not in all years (trend only for rows which are in both data sets)
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




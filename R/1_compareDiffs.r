### Hilfsfunktion zur Berechnung von cross-level group differences (fuer jk2.mean und jk2.table)
computeCrossLevel <- function ( jk2, cols, grpv, fun, cl_diffs, comp_type = NULL, tv = NULL) {
       ori <- jk2
    ### Achtung: rows nur auswaehlen, wenn imput aus jk2.mean stammt, bei jk2.table immer alles nehmen!
       if ( fun == "mean" ) {
            jk2  <- jk2[which(jk2[,"parameter"] %in% c("mean", "sd")),]
       }
       if ( fun %in% c("glm", "table") ) {
            jk2  <- jk2[which(!jk2[,"parameter"] %in% c("Nvalid", "Ncases", "R2", "R2nagel")),]
       }
    ### wenn die Funktions fuer TRends aufgerufen wird, gibt es nun mehrere trends (statt wie frueher bei 2 zeitpunkten nur einen)
    ### daher wird hier nach trendgruppen gesplittet, aber nur, wenn es trends gibt!! sonst wird eine Konstante zum Schleifenaufruf erzeugt
       if (is.null(tv)) {jk2[,"gruppierungsvariable"] <- "all"; tv <- "gruppierungsvariable"; isNullTv <- TRUE} else { isNullTv <- FALSE}
    ### cross-level diffs werden fuer die Gruppen 'comparison = NA' und 'comparison != NA' separat bestimmt, um
    ### die menge zu reduzieren und weil man die anderen nicht braucht
       retA<- by ( data = jk2, INDICES = jk2[,tv], FUN = function ( trenddat) {
              gr  <- data.frame ( a = trenddat[,"comparison"], b = is.na(trenddat[,"comparison"]))
              gr  <- paste(gr[,"a"], gr[,"b"], sep="_")
              ret <- do.call("rbind", by ( data = trenddat, INDICES = gr, FUN = function ( d ) {
    ### Hierarchy levels bestimmen
    ### Logik: alle mit einem nicht fehlenden Wert auf der Gruppenvariablen werden gegen die mit zwei fehlenden Werten verglichen
    ### alle mit keinem fehlenden Werte werden gegen die mit einem fehlenden Wert verglichen
    ### allgemein: alle mit x fehlenden Werten werden gegen alle mit x-1 fehlenden Werten verglichen
    ### ah == 'Anzahl hierarchien' == 'Anzahl Gruppenvariablen'
                     grp_log   <- data.frame ( !is.na(d[,grpv, drop=FALSE]) )
                     grp_log[,"sum"]    <- as.numeric(rowSums(grp_log))
                     d[,"sum"] <- grp_log[,"sum"]
    ### Abbruchkriterium definieren: keine Vergleiche berechnen, wenn nicht mindestens zwei unique elemente in 'd[,"sum"]'
                     if ( length(unique(d[,"sum"]))<=1) {
                          return(NULL)
                     }
    ### wenn cross-level diffs von group.diffs bestimmt werden sollen (comparison != NA), dann muessen die cl_diffs angepasst werden
                     if ( !is.na(d[1,"comparison"])) {
                          cl_diffs  <- as.list(data.frame ( combinat::combn(unique(d[,"sum"]),2)))
                     }
    ### loop over hierarchy levels ... Vergleichsrichtung wie in eatRep-Funktion festgelegt
                     ret <- do.call("rbind", lapply ( cl_diffs, FUN = function ( comp_vec ) {
                            redDF_1 <- d[which(d[, "sum"] %in% comp_vec),]
                            fac     <- by(data = redDF_1, INDICES = redDF_1[, "sum"], FUN = function ( x ) { unique(x[, "group"])})
    ### Achtung: nur ineinander geschachtelte CrossLevel-Differenzen werden gebildet
                            all_grp_list <- lapply(fac[[1]], function(high_lvl) {
    ### all comparisons if wholeGroup higher level
                                       hl_levels <- unlist(strsplit(high_lvl, "_"))
                                       if(identical(high_lvl, "wholeGroup") | substr(high_lvl, 1,9) == "all.group") {hl_levels <- "."}
    ### extrahiere alle in hoehere Ebene geschachtelten lower levels
                 # falls diese doch gebraucht werden:
                 # if(allCrossLvlDiffs) name_part <- "."
                                       all_ll_levels <- lapply(fac[[2]], FUN = function ( x) {unlist(strsplit(x, "_"))})
                                       ## Suche matches
                                       matching_levels <- sapply(all_ll_levels, FUN = function ( a ) {all(hl_levels %in% a)})
                                       low_lvl <- fac[[2]][matching_levels]
                                       if ( length(low_lvl)==0) {
    ### workaround: wenn cross-level diffs von group.diffs bestimmt werden sollen (comparison != NA), dann muessen 'low_lvl' anders gefunden werden
                                            low_lvl <- grep(reshape2::colsplit(string = paste0(hl_levels, collapse = "_"), pattern="___", names=c("a", "b"))[,"a"], fac[[2]], value = TRUE)
                                       }
                                       vgl     <- expand.grid(high_lvl, low_lvl)
    ### loop over comparison to be made!
                                       grp <- do.call("rbind", plyr::alply(as.matrix(vgl), .margins = 1, .fun = function(single_comp) {
                                              redDF_2 <- redDF_1[which(redDF_1[,"group"] %in% single_comp),]
    ### sort corresponding to order in cross.diff object (first number, second number), then that way the difference is calc
                                              redDF_2 <- redDF_2[c(which(redDF_2[, "sum"] == comp_vec[1]), which(redDF_2[, "sum"] == comp_vec[2])), ]
    ### calculate difference between all parameters
                                              newRows <- compareParameters(df_allP = redDF_2, grpv = grpv, fun = fun, comp_type = comp_type)
    ### Output aufbereiten fuer cross-level diffs von group.diffs
                                              return(newRows)
                                       }))
                                       return(grp)
                            })
                            all_grp <- do.call("rbind", all_grp_list)
                            return(all_grp)})
                     )
                     return(ret)
              })) })                                                            ### untre Zeile: Hotfix: dummyvariable loeschen falls es keine trends gibt
       if (isNullTv) {retA[[1]] <- retA[[1]][,-match("gruppierungsvariable", colnames(retA[[1]]))]}
       ret <- rbind(ori, do.call("rbind", retA))
       return(ret) }

########### Compare all parameters of a comparison (longformat)
### common function for crosslevel and trenddifferences
# comp_type only for naming in comparison variable!
compareParameters <- function(df_allP, grpv, fun, comp_type = NULL) {
  # SW: output ggf. sortieren
  df_allP<- df_allP[with(df_allP, order(parameter, group)),]
  # drop all irrelevant coefficients
  df_allP   <- df_allP[which(df_allP[,"coefficient"] %in% c("est", "se")),]
  ## for means: extract SDs
  if(identical(fun, "mean")) { df_sd <- df_allP[df_allP[, "parameter"] == "sd", ] }
  out <- by ( data = df_allP, INDICES = df_allP[,"parameter"], FUN = function ( df ) {
    stopifnot(nrow(df) == 4 || nrow(df)==2)                                     ### checks
    #if(df$group[1] == "female") browser()
    df <- df[order(df$sum, decreasing = FALSE), ]                                ### Direction of crossDiff: Higher vs Lower (eg country vs all)
    mea <- diff(df[which(df[,"coefficient"] == "est"),"value"])                 ### compute mean difference
    if ( length(mea) ==0 || is.na(mea)) {return(NULL)}
	if ( !"se" %in% df[, "coefficient"]) {
         cat(paste0( "   Warning: No standard error for parameter '",unique(df[,"parameter"]),"'. Cannot compute standard errors and p value for difference between '",df[1,"group"],"' and '",df[2,"group"],"'.\n"))
         se <- pval <- NA
    }  else  {
         se  <- sqrt(sum(df[which(df[,"coefficient"] == "se"),"value"]^2))      ### compute SE
         pval<- 2*pnorm(abs(mea/se), lower.tail=FALSE)
    }
    es  <- NA
    if (  fun == "mean" && df[1,"parameter"] == "mean" && nrow(df_sd)>0) {
      sd_wide <- reshape2::dcast(df_sd[which(df_sd[,"coefficient"] == "est"),], group~parameter, value.var = "value")
   ### achtung: wenn cross-differences fuer adjusted means gemacht werden, gibt es manchmal keine Standardabweichung fuer Gruppenmittelwerte
   ### sd_wide hat dann nur eine Zeile. wenn das so ist, kann keine effektstaerke berechnet werden   
      if ( nrow(sd_wide) > 1) { es  <- mea / sqrt(0.5*sum(sd_wide[,"sd"]^2)) }
    }
    if(nrow(df) == 2 ) {                                                        ### uebler Hotfix: ohne Jackknife gibt es keine Se fuer SDs
       ret <- rbind (df, df)                                                    ### 'df' hat dann nur 2 statt 4 Zeilen
    }  else  {
       ret <- df
    }
    # if no further specified use old comparison type (for trends)
    if(is.null(comp_type)) {
      ret[,"comparison"] <-  df[1, "comparison"]
    } else  {
      ret[,"comparison"] <- comp_type
    }
    ret[,"coefficient"] <- c("est","se", "p", "es")
    ret[,"value"]       <- c(mea, se, pval, es)

    ### group column: consistent ordering of higher hierarchy level
    hierarchy_levels <- unique(df[order(df$sum, decreasing = TRUE), "group"])
    ret[,"group"]       <- paste(hierarchy_levels, collapse=".vs.")
    ret[,"sum"]         <- NULL
    # set group variable to NA for group that is compared against complete group (other remain)
    for(i in grpv) {
      if(any(is.na(ret[, i]))) ret[, i] <- NA
    }
    # drop empty rows
    ret <- ret[which(!is.na(ret[,"value"])),]
    ret
  })
  return(do.call("rbind", out))
}
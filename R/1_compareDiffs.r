### Hilfsfunktion zur Berechnung von cross-level group differences (fuer jk2.mean und jk2.table)
computeCrossLevel <- function ( jk2, cols, grpv, fun, cl_diffs, comp_type = NULL) {
       ori <- jk2
    ### Achtung: rows nur auswaehlen, wenn imput aus jk2.mean stammt, bei jk2.table immer alles nehmen!
       if ( fun == "mean" ) {
            jk2  <- jk2[which(jk2[,"parameter"] %in% c("mean", "sd")),]
       }
       if ( fun == "glm" ) {
            jk2  <- jk2[which(!jk2[,"parameter"] %in% c("Nvalid", "Ncases", "R2", "R2nagel")),]
       }
    ### cross-level diffs werden fuer die Gruppen 'comparison = NA' und 'comparison != NA' separat bestimmt, um
    ### die menge zu reduzieren und weil man die anderen nicht braucht
       gr  <- data.frame ( a = jk2[,"comparison"], b = is.na(jk2[,"comparison"]))
       gr  <- paste(gr[,"a"], gr[,"b"], sep="_")
       ret <- do.call("rbind", by ( data = jk2, INDICES = gr, FUN = function ( d ) {
    ### Hierarchy levels bestimmen
    ### Logik: alle mit einem nicht fehlenden Wert auf der Gruppenvariablen werden gegen die mit zwei fehlenden Werten vergleichen
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
                   cl_diffs  <- as.list(data.frame ( combn(unique(d[,"sum"]),2)))
              }
    ### loop over hierarchy levels ... Vergleichsrichtung wie in eatRep-Funktion festgelegt
              ret <- do.call("rbind", lapply ( cl_diffs, FUN = function ( comp_vec ) {
                     redDF_1 <- d[which(d[, "sum"] %in% comp_vec),]
                     fac     <- by(data = redDF_1, INDICES = redDF_1[, "sum"], FUN = function ( x ) { unique(x[, "group"])})
    ### Achtung: nur ineinander geschachtelte CrossLevel-Differenzen werden gebildet
                     all_grp <- do.call("rbind", lapply(fac[[1]], function(high_lvl) {
    ### all comparisons if wholeGroup higher level
                                name_part <- high_lvl
                                if(identical(high_lvl, "wholeGroup") | substr(high_lvl, 1,9) == "all.group") {name_part <- "."}
    ### extrahiere alle in hoehere Ebene geschachtelten lower levels
          # falls diese doch gebraucht werden:
          # if(allCrossLvlDiffs) name_part <- "."
                                low_lvl <- grep(name_part, fac[[2]], value = TRUE)
                                if ( length(low_lvl)==0) {
    ### workaround: wenn cross-level diffs von group.diffs bestimmt werden sollen (comparison != NA), dann muessen 'low_lvl' anders gefunden werden
                                     low_lvl <- grep(colsplit(string = name_part, pattern="___", names=c("a", "b"))[,"a"], fac[[2]], value = TRUE)
                                }
                                vgl     <- expand.grid(high_lvl, low_lvl)
    ### loop over comparison to be made!
                                grp <- do.call("rbind", alply(as.matrix(vgl), .margins = 1, .fun = function(single_comp) {
                                       redDF_2 <- redDF_1[which(redDF_1[,"group"] %in% single_comp),]
    ### sort corresponding to order in cross.diff object (first number, second number), then that way the difference is calc
                                       redDF_2 <- redDF_2[c(which(redDF_2[, "sum"] == comp_vec[1]), which(redDF_2[, "sum"] == comp_vec[2])), ]
    ### calculate difference between all parameters
                                       newRows <- compareParameters(df_allP = redDF_2, grpv = grpv, fun = fun, comp_type = comp_type)
    ### Output aufbereiten fuer cross-level diffs von group.diffs
                                       return(newRows)
                                }))
                                return(grp)
                     }))
                     return(all_grp)
              }))
              return(ret)
       }))
       ret <- rbind(ori, ret)
       return(ret) }

########### Compare all parameters of a comparison (longformat)
### common function for crosslevel and trenddifferences
# comp_type only for naming in comparison variable!
compareParameters <- function(df_allP, grpv, fun, comp_type = NULL) {
  # drop all irrelevant coefficients
  df_allP   <- df_allP[which(df_allP[,"coefficient"] %in% c("est", "se")),]
  ## for means: extract SDs
  if(identical(fun, "mean")) { df_sd <- df_allP[df_allP[, "parameter"] == "sd", ] }
  out <- by ( data = df_allP, INDICES = df_allP[,"parameter"], FUN = function ( df ) {
    stopifnot(nrow(df) == 4 || nrow(df)==2)                                     ### checks
    mea <- diff(df[which(df[,"coefficient"] == "est"),"value"])                 ### compute mean difference
    if ( !"se" %in% df[, "coefficient"]) {
         cat(paste0( "   Warning: No standard error for parameter '",unique(df[,"parameter"]),"'. Cannot compute standard errors and p value for difference between '",df[1,"group"],"' and '",df[2,"group"],"'.\n"))
         se <- pval <- NA
    }  else  {
         se  <- sqrt(sum(df[which(df[,"coefficient"] == "se"),"value"]^2))      ### compute SE
         pval<- 2*pnorm(abs(mea/se), lower=FALSE)
    }
    es  <- NA
    if (  fun == "mean" && df[1,"parameter"] == "mean" && nrow(df_sd)>0) {
      sd_wide <- dcast(df_sd[which(df_sd[,"coefficient"] == "est"),], group~parameter, value.var = "value")
      es  <- mea / sqrt(0.5*sum(sd_wide[,"sd"]^2))
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
    ret[,"group"]       <- paste(unique(df[,"group"]), collapse=".vs.")
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

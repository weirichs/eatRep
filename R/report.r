report <- function ( jk2.out, trendDiffs = FALSE, add=list(), exclude = c("Ncases", "NcasesValid", "var"), printGlm = FALSE, round = TRUE, digits = 3, printDeviance = FALSE) {
          if(is.null(jk2.out)) {return(NULL)}
    ### vorab: alte 'dG'-Funktion zum Anzeigen der Regressionsergebnisse implementieren
          if ( length(grep("glm", as.character(jk2.out[["resT"]][[1]][1,"modus"]))) ==1 ) {
               if ( printGlm == TRUE ) { dG(jk2.out, digits = digits, printDeviance = printDeviance ) }
          }
    ### 1. Input extrahieren: diese Variablen dann spaeter an Einzelfunktionen weitergeben!
          jk2      <- jk2.out[["resT"]]
          tv       <- jk2.out[["allNam"]][["trend"]]
          cols     <- c("group", "depVar",  "modus", "parameter")
          grpv     <- setdiff(setdiff(colnames(jk2[[1]]), cols), c("comparison", "coefficient", "value", tv))
          grp_by   <- jk2.out[["allNam"]][["group.differences.by"]]
          cl_diffs <- jk2.out[["allNam"]][["cross.differences"]]
          funs     <- c("mean", "table", "quantile", "glm")
          fun      <- funs [ which( unlist(lapply(funs, FUN = function ( f ) { length(grep(f, jk2[[1]][1,"modus"]))})) > 0) ]
    ### 2. cross-level diffs bestimmen: ueberschreibt bzw. erweitert das Objekt 'jk2' ... Achtung: sind nur fuer "mean" oder "table" erlaubt
          if ( is.list(cl_diffs) ) {
               jk2 <- lapply(jk2, FUN = function (df) {computeCrossLevel (df, cols=cols, grpv = grpv, fun = fun, cl_diffs = cl_diffs, comp_type = "crossDiff")})
               jk2 <- lapply(jk2, FUN = function (df) {                         ### in spalte "comparison" 'crossDiff_of_groupDiff' eintragen, falls in Spalte "group" 3x ".vs." steht
                      spl <- strsplit(df[,"group"], ".vs.")
                      ind <- which(sapply(spl, length)==4)
                      if ( length(ind)>0) {
                           df[ind,"comparison"] <- "crossDiff_of_groupDiff"     ### das mittlere ".vs." gross schreiben
                           df[ind,"group"]      <- unlist(lapply(spl[ind], FUN = function ( z ) {paste(z[1], ".vs.", z[2], ".VS.", z[3], ".vs.", z[4], sep="")}))
                      }
                      return(df)})
          }
    
    ### 2.b) SE correction durchfuehren (siehe Paper Weirich & Hecht)
          if(!is.null(jk2.out[["SE_correction"]]) && !is.null(jk2.out[["SE_correction"]][[1]])) {
            ## checks, ob Vergleiche dabei, fuer die keine Korrektur verfuehgbar ist
            if(length(which(jk2[[1]][["comparison"]] == "crossDiff_of_groupDiff")) > 0 ) {
              warning("Standard error correction for 'crossDiff_of_groupDiff' is currently not supported.")
            }
            mult_hierarchy <- any(unlist(lapply(jk2.out$allNam$cross.differences, function(x) x[2] - x[1] != 1)))
            if(mult_hierarchy) warning("Standard error correction for crossDifferences across multiple hierarchy levels is currently not supported.")
            
            # browser()
            ## correction durchfuehren
            jk2 <- lapply(jk2, function(jk2_single) {
              seCorrect(SE_correction = jk2.out[["SE_correction"]], jk2 = jk2_single, grpv = grpv)
            })
          }
          
    ### 3. Trend bestimmen
          if ( !is.null(tv) ) {
               jk2 <- computeTrend(jk2 = jk2, le = jk2.out[["le"]], tv = tv, fun = fun)
          } else {
               jk2 <- jk2[[1]]
          }
    ### 4. Trend-Differences bestimmen ... Achtung: sind nur fuer "mean" oder "table" erlaubt
          if ( !is.null(tv) && trendDiffs ) {
               jk2 <- computeTrendDiffs(jk2 = jk2, grpv = grpv, tv = tv, grp_by = grp_by, fun = fun, cl_diffs = cl_diffs)
          }
    ### 5. 'add' ergaenzen, falls gewuenscht
          if ( length(add)>0) {
               if(!all(nchar(names(add))>0)) { stop("'add' must be named.")}    ### necessary checks
               if(length(names(add)) != length(unique(names(add)))) { stop("Duplicated names of 'add' are not allowed.")}
               if(!all(sapply(add, length) == 1)) {stop("All elements of 'add' must be of length 1.")}
               if(!all(sapply(add, class) == "character")) {stop("All elements of 'add' must be of class 'character'.")}
               dopp<- names(add) %in% colnames(jk2)
               ind <- which(dopp==TRUE)
               if ( length( ind ) > 0 ) {stop(paste0("Following names of 'add' are not allowed: '",paste(names(add)[ind], collapse = "', '"), "'."))}
               for ( u in names(add)) {jk2[,u] <- add[[u]]}
          }
    ### 6. reshapen
          spltVar  <- c("coefficient", tv)                                      ### split-Variable
          if ( length(exclude)>0) {
               weg <- which(jk2[,"parameter"] %in% exclude)
               if ( length(weg)>0) {
                    jk2 <- jk2[-weg,]
               }
          }                                                                     ### was muss in die Spalten? das haengt davon ab, ob es einen Trend gibt
          frml     <- as.formula(paste0("... ~ ", paste(spltVar,collapse=" + ") ) )
          jk2wide  <- dcast(data = jk2, formula = frml, value.var = "value")
    ### runden, falls gewuenscht
          if ( round == TRUE) {
               coln<- which(sapply(jk2wide, is.numeric))
               if ( length(coln)>0) {
                    for ( i in coln) { jk2wide[,i] <- round(jk2wide[,i], digits = digits)}
               }
          }
          return(jk2wide)}

addSig <- function ( dat , groupCols = NULL , allNam = NULL ) {
          if(is.null(groupCols)) {groupCols <- c("group", "parameter")}
          dat <- do.call("rbind", by ( data = dat, INDICES = dat[,groupCols], FUN = function ( x ) {
                 z  <- x[which(x[,"coefficient"] %in% c("est", "se")),]
                 if ( nrow(z) > 2) {cat("Fehler. x muss maximal 2 zeilen haben.\n")}
                 if ( nrow(z) == 2 ) {
                      y  <- z[1,]                                               ### dazu relevante spalten identifizieren, nach denen gesplittet werden muss
                      y[["coefficient"]] <- "p"                                 ### erste Zeile von x duplizieren und relevante Werte ersetzen
                      y[["value"]]       <- 2*pnorm(abs(z[which(z[,"coefficient"] == "est"),"value"] / z[which(z[,"coefficient"] == "se"),"value"]), lower=FALSE)
                      x  <- rbind ( x, y)                                       ### Achtung: Signifikanzwert wird hier noch nach numerisch transformiert, muss zurueckgewandelt werden
                 }
                 return(x)}))                                                   ### untere Zeile: wenn 'table' ueber 'jk2.mean' gewrappt wurde, muessen hier die parameterbezeichnungen geaendert werden
          return(dat)}


seCorrect <- function( SE_correction, jk2, grpv ) {
  # stop("SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
  UseMethod("seCorrect")
}

seCorrect.old <- function( SE_correction, jk2, grpv ) {
  jk2
}

## an der falschen Stelle! muss vor Trends passieren!
seCorrect.wec_se_correction <- function( SE_correction, jk2, grpv ) {
  #stop("SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
  
  ## Trend or no trend
  # if(is.null(jk2[["year"]])) browser()
  year <- as.character(unique(jk2[["year"]]))
  if(length(year) == 0) year <- 1

  ## dissect results so far
  no_cross_diff <- jk2[is.na(jk2$comparison) | jk2$comparison != "crossDiff", ]
  cross_diff <- jk2[which(jk2$comparison == "crossDiff"), ]
  
  ## mehr Fragen
  # _ in Faktorlevels erlaubt?
  
  ## Fragen
  # crossDiff of groupDiff? -> wird hinten angestellt! vlt. warning/cat, dass es dafuer nicht implementiert ist
  
  # Vgl. ueber mehrere lvls (zB Male in Berlin vs. Deutschland): da hat Sebastian keine Lust das zu implementieren,
  # findet er inhaltlich keine spannende Fragestellung!
  # vlt. Warnmeldung?
  
  # generell nochmal ueberlegen: einfacher zu ersetzn oder in crossDiff funktion einzubauen?
  
  # if(length(SE_correction) > 2) browser()
  for(i in seq_along(SE_correction)) {
    output <- SE_correction[[i]][["resT"]][[year]]
    
    rows <- length(SE_correction[[i]][["vgl"]][["groups.divided.by"]])
    single_grpv <- as.character(SE_correction[[i]][["vgl"]][["groups.divided.by"]])[rows]
    
    #if(grepl(" \\+ ", single_grpv)) browser()
    if(grepl(" \\+ ", single_grpv)) {
      single_grpv <- gsub(SE_correction[[i]][["vgl"]][["groups.divided.by"]][rows-1], "", single_grpv)
      single_grpv <- gsub(" \\+ ", "", single_grpv)
    }
    
    ### Achtung: ueberpruefen, ob allNam$independent nicht die Info einfacher bereit haelt
    
    ## huihui, Achtung: jenachdem ob Trend oder nicht, heißen die Spalten ja anders!
    # also eher nicht die reporting Funktion verwenden
    SEs <- output[!output$parameter %in% c("(Intercept)", "Nvalid", "R2") & output$coefficient %in% c("se", "p"), c("parameter", "value", "coefficient")]
    SEs[, "parameter"] <- gsub(single_grpv, "", SEs[, "parameter"])
    SEs <- as.data.frame(tidyr::pivot_wider(SEs, names_from = "coefficient", values_from = "value"))
    
    for(param in SEs[["parameter"]]) {
      ## if reference level is the whole group
      if(identical(SE_correction[[i]]$refGrp, "all")) {
        grp_regexp <- paste0("^", param, "\\.vs")
        # funktioniert, wenn die Gesamtheit die Vergleichsgruppe ist
        cross_diff[cross_diff$parameter == "mean" & grepl(grp_regexp, cross_diff$group) & cross_diff$coefficient == "se", 
                   "value"] <- SEs[SEs[, "parameter"] == param, "se"]
        cross_diff[cross_diff$parameter == "mean" & grepl(grp_regexp, cross_diff$group) & cross_diff$coefficient == "p", 
                   "value"] <- SEs[SEs[, "parameter"] == param, "p"]
      } else { ## if reference level is a subgroup
        grp_regexp <- paste0("_", param, "$")
        
        col_name <- SE_correction[[i]]$refGrp[, "groupName"]
        col_level <- SE_correction[[i]]$refGrp[, "groupValue"]
        # if(col_level == "LandA") browser()
        
        cross_diff[which(cross_diff[[col_name]] == col_level & grepl(grp_regexp, cross_diff$group) & cross_diff$coefficient == "se"), 
                   "value"] <- SEs[SEs[, "parameter"] == param, "se"]
        cross_diff[which(cross_diff[[col_name]] == col_level & grepl(grp_regexp, cross_diff$group) & cross_diff$coefficient == "p"), 
                   "value"] <- SEs[SEs[, "parameter"] == param, "p"]
      }
    }
  }
  
  ### next: trend, groupdiffs checken, ob hier SEs korrekt beruecksichtigt werden
  
  rbind(no_cross_diff, cross_diff)
}



seCorrect.pisa_se_correction <- function( SE_correction, jk2, grpv ) {
  stop("SE correction has not been implemented yet. Use crossDiffSE = 'old'.")
}



report <- function ( jk2.out, trendDiffs = FALSE, add=list(), exclude = c("Ncases", "NcasesValid", "var"), printGlm = FALSE, digits = 3, printDeviance = FALSE) {
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
    ### to do: gesamtfiltervariable bestimmen (das was Felix/Steffi haben wollten)
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

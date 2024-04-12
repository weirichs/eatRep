checkLEs <- function(filePaths, leDF) {
  checkmate::assert_file_exists(filePaths)
  checkmate::assert_character(filePaths, min.len = 2)
  leDF <- eatTools::makeDataFrame(leDF)
    ### erstmal nur das linking error objekt an sich checken: check variable names
       allV  <- list(trendLevel1 = "trendLevel1", trendLevel2 = "trendLevel2", parameter = "parameter", linkingError="linkingError", depVar = "depVar")
       allN  <- lapply(allV, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = leDF, variable=ii, warnIfMissing = TRUE)})
    ### Namen in GADSdat data bases
       namlis<- lapply(filePaths, eatGADS::namesGADS)
    ### all Linking error variables in trend gads data bases (without LE_)?
       depVar<- unique(leDF[,allN[["depVar"]]])
       dep_notIn_nam_list <- lapply(seq_along(namlis), function(i) {
           nam  <- namlis[[i]]
           dep_notIn_nam <- setdiff(depVar, unlist(nam))
           if(length(dep_notIn_nam) > 0) {message("The following variables have linking errors but are not variables in data base ",  i, ": '", paste0(dep_notIn_nam, collapse = "', '"),"'")}
           return(dep_notIn_nam)  })
    ### Anzahl der Trendzeitpunkte (trend levels)
       tLevel<- unique(unlist(leDF[,c(allN[["trendLevel1"]], allN[["trendLevel2"]])]))
       if ( length(tLevel) != length(filePaths)) {warning(paste0("Number of trend levels do not match: Expect ",length(filePaths)," trend levels in data base ('filePaths' has length ",length(filePaths),"). ",length(tLevel), " trend levels for linking errors ('",paste(tLevel, collapse="', '"),"') found."))}
       return(list(dep_notIn_nam = dep_notIn_nam_list)) }
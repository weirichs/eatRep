checkLEs <- function(filePaths, leDF) {
  nam_list <- lapply(filePaths, eatGADS::namesGADS)
  namLE <- names(leDF)
  
  ## check variable names?
  if(!identical(namLE, c("trendLevel1", "trendLevel2", "depVar", "domain", "parameter", "linkingError"))) stop("Incorrect variable names in 'leDF'.")
  
  ## all Linking error variables in trend gads data bases (without LE_)?
  dep_variables <- leDF$depVar
  dep_notIn_nam_list <- lapply(seq_along(nam_list), function(i) {
    nam <- nam_list[[i]]
    dep_notIn_nam <- setdiff(dep_variables, unlist(nam))
    if(length(dep_notIn_nam) > 0) message("The following variables have linking errors but are not variables in data base ",  i, ": ",
                                          paste0(dep_notIn_nam, collapse = ", "))
    dep_notIn_nam
  })
  
  
    ## potenziell: Domain checken? Aber dafuer muesste man Variablenname festlegen
  # domain_variables <- leDF$domain
  
  list(dep_notIn_nam = dep_notIn_nam_list)
}
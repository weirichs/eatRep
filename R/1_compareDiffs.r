computeCrossLevel <- function ( jk2, cols, grpv, fun, cl_diffs, comp_type = NULL) {
       ori <- jk2
       if ( fun == "mean" ) {
            jk2  <- jk2[which(jk2[,"parameter"] %in% c("mean", "sd")),]
       }
       if ( fun == "glm" ) {
            jk2  <- jk2[which(!jk2[,"parameter"] %in% c("Nvalid", "Ncases", "R2", "R2nagel")),]
       }
       gr  <- data.frame ( a = jk2[,"comparison"], b = is.na(jk2[,"comparison"]))
       gr  <- paste(gr[,"a"], gr[,"b"], sep="_")
       ret <- do.call("rbind", by ( data = jk2, INDICES = gr, FUN = function ( d ) {
              grp_log   <- data.frame ( !is.na(d[,grpv, drop=FALSE]) )
              grp_log[,"sum"]    <- as.numeric(rowSums(grp_log))
              d[,"sum"] <- grp_log[,"sum"]
              if ( length(unique(d[,"sum"]))<=1) {
                   return(NULL)
              }
              if ( !is.na(d[1,"comparison"])) {
                   cl_diffs  <- as.list(data.frame ( combinat::combn(unique(d[,"sum"]),2)))
              }
              ret <- do.call("rbind", lapply ( cl_diffs, FUN = function ( comp_vec ) {
                     redDF_1 <- d[which(d[, "sum"] %in% comp_vec),]
                     fac     <- by(data = redDF_1, INDICES = redDF_1[, "sum"], FUN = function ( x ) { unique(x[, "group"])})
                     all_grp_list <- lapply(fac[[1]], function(high_lvl) {
                                hl_levels <- unlist(strsplit(high_lvl, "_"))
                                if(identical(high_lvl, "wholeGroup") | substr(high_lvl, 1,9) == "all.group") {hl_levels <- "."}
                                all_ll_levels <- lapply(fac[[2]], FUN = function ( x) {unlist(strsplit(x, "_"))})
                                matching_levels <- sapply(all_ll_levels, FUN = function ( a ) {all(hl_levels %in% a)})
                                low_lvl <- fac[[2]][matching_levels]
                                if ( length(low_lvl)==0) {
                                     low_lvl <- grep(reshape2::colsplit(string = paste0(hl_levels, collapse = "_"), pattern="___", names=c("a", "b"))[,"a"], fac[[2]], value = TRUE)
                                }
                                vgl     <- expand.grid(high_lvl, low_lvl)
                                grp <- do.call("rbind", plyr::alply(as.matrix(vgl), .margins = 1, .fun = function(single_comp) {
                                       redDF_2 <- redDF_1[which(redDF_1[,"group"] %in% single_comp),]
                                       redDF_2 <- redDF_2[c(which(redDF_2[, "sum"] == comp_vec[1]), which(redDF_2[, "sum"] == comp_vec[2])), ]
                                       newRows <- compareParameters(df_allP = redDF_2, grpv = grpv, fun = fun, comp_type = comp_type)
                                       return(newRows)
                                }))
                                return(grp)
                     })
                     all_grp <- do.call("rbind", all_grp_list)
                     return(all_grp)})
              )
              return(ret)
       }))
       ret <- rbind(ori, ret)
       return(ret) }

compareParameters <- function(df_allP, grpv, fun, comp_type = NULL) {
  df_allP<- df_allP[with(df_allP, order(parameter, group)),]
  df_allP   <- df_allP[which(df_allP[,"coefficient"] %in% c("est", "se")),]
  if(identical(fun, "mean")) { df_sd <- df_allP[df_allP[, "parameter"] == "sd", ] }
  out <- by ( data = df_allP, INDICES = df_allP[,"parameter"], FUN = function ( df ) {
    stopifnot(nrow(df) == 4 || nrow(df)==2)                                     
    df <- df[order(df$sum, decreasing = FALSE), ]                                
    mea <- diff(df[which(df[,"coefficient"] == "est"),"value"])                 
    if ( length(mea) ==0 || is.na(mea)) {return(NULL)}
	if ( !"se" %in% df[, "coefficient"]) {
         cat(paste0( "   Warning: No standard error for parameter '",unique(df[,"parameter"]),"'. Cannot compute standard errors and p value for difference between '",df[1,"group"],"' and '",df[2,"group"],"'.\n"))
         se <- pval <- NA
    }  else  {
         se  <- sqrt(sum(df[which(df[,"coefficient"] == "se"),"value"]^2))      
         pval<- 2*pnorm(abs(mea/se), lower.tail=FALSE)
    }
    es  <- NA
    if (  fun == "mean" && df[1,"parameter"] == "mean" && nrow(df_sd)>0) {
      sd_wide <- reshape2::dcast(df_sd[which(df_sd[,"coefficient"] == "est"),], group~parameter, value.var = "value")
      if ( nrow(sd_wide) > 1) { es  <- mea / sqrt(0.5*sum(sd_wide[,"sd"]^2)) }
    }
    if(nrow(df) == 2 ) {                                                        
       ret <- rbind (df, df)                                                    
    }  else  {
       ret <- df
    }
    if(is.null(comp_type)) {
      ret[,"comparison"] <-  df[1, "comparison"]
    } else  {
      ret[,"comparison"] <- comp_type
    }
    ret[,"coefficient"] <- c("est","se", "p", "es")
    ret[,"value"]       <- c(mea, se, pval, es)

    hierarchy_levels <- unique(df[order(df$sum, decreasing = TRUE), "group"])
    ret[,"group"]       <- paste(hierarchy_levels, collapse=".vs.")
    ret[,"sum"]         <- NULL
    for(i in grpv) {
      if(any(is.na(ret[, i]))) ret[, i] <- NA
    }
    ret <- ret[which(!is.na(ret[,"value"])),]
    ret
  })
  return(do.call("rbind", out))
}

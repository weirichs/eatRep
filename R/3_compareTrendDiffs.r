computeTrendDiffs <- function(jk2, grpv, tv, grp_by, fun,  cl_diffs) {
  ori <- jk2
  jk2_trend <- jk2[jk2[tv] == "trend", ]
  jk2_crossDiff <- NULL
  if(is.list(cl_diffs)) {
    jk2_crossDiff <- computeCrossLevel(jk2 = jk2_trend, grpv = grpv, fun = fun,  cl_diffs = cl_diffs,
                                  comp_type = "trendDiff_cross")
  }
  jk2_groupDiff <- NULL
  if(!is.null(grp_by)) jk2_groupDiff <- trendGroupDiffs(jk2 = jk2_trend, grpv = grpv, grp_by = grp_by, fun = fun)
  jk2_all <- rbind(ori, jk2_crossDiff, jk2_groupDiff)
  return(unique(jk2_all))
}

trendGroupDiffs <- function(jk2, grpv, grp_by, fun) {
  ori <- jk2
  jk2_groups <- jk2[jk2[, "group"] != "wholeGroup", ]
  jk2_groups <- jk2_groups[!is.na(jk2_groups[, grp_by]), ]
  jk2 <- jk2_groups
  if ( fun == "mean" ) {
    rows <- intersect(which(jk2[,"parameter"] %in% c("mean", "sd")), which(is.na(jk2[,"comparison"])))
    jk2 <- jk2[rows,]
  }  else  {
    jk2 <- jk2[which(is.na(jk2[,"comparison"])),]
  }
  grp_log <- data.frame ( !is.na(jk2[,grpv, drop=FALSE]) )
  grp_log[,"sum"]    <- as.numeric(rowSums(grp_log))
  stopifnot(!"sum" %in% colnames(jk2))
  jk2[,"sum"] <- grp_log[,"sum"]
  grp_notby <- grpv[!grpv  %in% grp_by]
  ret <- by(jk2, INDICES = jk2[, "sum"], FUN = function(redDF_1) {
    fac_lvls <- redDF_1[, grp_notby]
    fac_lvls[is.na(fac_lvls)] <- "none"
    if(length(fac_lvls) == 0) { fac_lvls <- rep("none", nrow(redDF_1)) }        
    grp_set <- by(redDF_1, INDICES = fac_lvls, function(redDF_2) {
      redDF_2 <- redDF_2[order(redDF_2[, "group"]), ]
      sorted_levels <- sort(unique(redDF_2[, "group"]))
      vgl <- t(combinat::combn(sorted_levels, 2))
      grp <- do.call("rbind", plyr::alply(as.matrix(vgl), .margins = 1, .fun = function(single_comp) {
        aus <- redDF_2[which(redDF_2[,"group"] %in% single_comp),]
        newRows <- compareParameters(df_allP = aus, grpv = grpv, fun = fun,
                                     comp_type = "trendDiff_group")
        fac_val <- vector(length = 0)
        count <- 1
        for(n in grp_notby) {
          if(!all(is.na(newRows[, n]))) {
            fac_val[count] <- paste(n, "=", unique(newRows[, n]), sep ="" )
            count <- count + 1
          }
        }
        if(length(fac_val) == 0 ) {
          fac_val <- "all.group=1"}
        else fac_val <- paste(fac_val, collapse = ", ")
        newRows[, grp_by] <- paste(unique(newRows[,grp_by]), collapse=".vs.")
        newRows[, "group"] <- paste(fac_val, "____", newRows[, grp_by], sep = "")
        return(newRows)
      }))
      return(grp)
    })
    grp_set <- do.call("rbind", grp_set)
    return(grp_set)
  })
  ret <- do.call(rbind, ret)
  return(rbind(ori, ret))
}

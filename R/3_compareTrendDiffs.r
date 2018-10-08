computeTrendDiffs <- function(jk2, grpv, tv, grp_by, fun,  cl_diffs) {
  ori <- jk2
  # pick only trends
  jk2_trend <- jk2[jk2[tv] == "trend", ]
  ### 1.: cross level diffs (only if possible)
  jk2_crossDiff <- NULL
  if(is.list(cl_diffs)) {
    jk2_crossDiff <- computeCrossLevel(jk2 = jk2_trend, grpv = grpv, fun = fun,  cl_diffs = cl_diffs,
                                  comp_type = "trendDiff_cross")
  }
  ### 2.: group level diffs (fuer group.diff.by variable)
  jk2_groupDiff <- NULL
  if(!is.null(grp_by)) jk2_groupDiff <- trendGroupDiffs(jk2 = jk2_trend, grpv = grpv, grp_by = grp_by, fun = fun)
  ### return results
  jk2_all <- rbind(ori, jk2_crossDiff, jk2_groupDiff)
  return(unique(jk2_all))
}

trendGroupDiffs <- function(jk2, grpv, grp_by, fun) {
  ori <- jk2
  # drop wholeGroup and rows that don't contain group difference variable
  jk2_groups <- jk2[jk2[, "group"] != "wholeGroup", ]
  jk2_groups <- jk2_groups[!is.na(jk2_groups[, grp_by]), ]
  jk2 <- jk2_groups
  ### Achtung: rows nur auswaehlen, wenn imput aus jk2.mean stammt, bei jk2.table immer alles nehmen!
  if ( fun == "mean" ) {
    rows <- intersect(which(jk2[,"parameter"] %in% c("mean", "sd")), which(is.na(jk2[,"comparison"])))
    jk2 <- jk2[rows,]
  }  else  {
    jk2 <- jk2[which(is.na(jk2[,"comparison"])),]
  }
  grp_log <- data.frame ( !is.na(jk2[,grpv, drop=FALSE]) )
  grp_log[,"sum"]    <- as.numeric(rowSums(grp_log))
  ### hier den Test auch raus???
  stopifnot(!"sum" %in% colnames(jk2))
  jk2[,"sum"] <- grp_log[,"sum"]
  grp_notby <- grpv[!grpv  %in% grp_by]
  # loop over hierarchy levels
  ret <- by(jk2, INDICES = jk2[, "sum"], FUN = function(redDF_1) {
    # for highest level: convert NA to none so by loops works over it normally
    fac_lvls <- redDF_1[, grp_notby]
    fac_lvls[is.na(fac_lvls)] <- "none"
    if(length(fac_lvls) == 0) { fac_lvls <- rep("none", nrow(redDF_1)) }        ### if no other group variables, initialize vector to by-loop over
    grp_set <- by(redDF_1, INDICES = fac_lvls, function(redDF_2) {
      # sort levels of group (in comparison object and data frame!) to stay consistent with jacknife.mean (Sebastian)
      redDF_2 <- redDF_2[order(redDF_2[, "group"]), ]
      sorted_levels <- sort(unique(redDF_2[, "group"]))
      # combine all combinations
      vgl <- t(combn(sorted_levels, 2))
      # loop over comparison to be made!
      grp <- do.call("rbind", alply(as.matrix(vgl), .margins = 1, .fun = function(single_comp) {
        aus <- redDF_2[which(redDF_2[,"group"] %in% single_comp),]
        # calculate difference between all parameters
        newRows <- compareParameters(df_allP = aus, grpv = grpv, fun = fun,
                                     comp_type = "trendDiff_group")
        # formatting group variables (make it comparable to groupDiff!)
        # like: sex=female____LandA.vs.LandB, sex=male, mig=0____LandA.vs.LandB
        fac_val <- vector(length = 0)
        count <- 1
        for(n in grp_notby) {
          if(!all(is.na(newRows[, n]))) {
            fac_val[count] <- paste(n, "=", unique(newRows[, n]), sep ="" )
            count <- count + 1
          }
        }
        # combine into 1 expression
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
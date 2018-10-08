computeTrend <- function(jk2, tv, le, fun) {
  # bind yearwise results
  jk2_all <- rbind(jk2[[1]], jk2[[2]])
  jk2_bind <- jk2_all
  # checks: die selben Zeilen in beiden Jahren?? (fuer GLMs insbesondere testen!)
  if(length(unique(jk2[[1]]$group)) != length(unique(jk2[[2]]$group)) ||
     sort(unique(jk2[[1]]$group)) != sort(unique(jk2[[2]]$group))) {
    uneven1 <- setdiff(unique(jk2[[1]]$group), unique(jk2[[2]]$group))
    uneven2 <- setdiff(unique(jk2[[2]]$group), unique(jk2[[1]]$group))
    if(length(uneven1) > 0) cat("Categories", uneven1, "are missing in year", jk2[[2]][1, tv] ,". \n")
    if(length(uneven2) > 0) cat("Categories", uneven2, "are missing in year", jk2[[1]][1, tv] ,". \n")
    # remove rows which are not in both years (trend only for rows which are in both data sets)
    jk2_bind <- jk2_all[!(jk2_all$group %in% c(uneven1, uneven2)), ]
  }
  # special for mean: select only mean and sd, reshape le
  if(identical(fun, "mean")) {
    jk2_bind <- jk2_bind[jk2_bind[["parameter"]] %in% c("mean", "sd"), ]
  }
  if ( fun == "glm" ) {
    jk2_bind <- jk2_bind[which(!jk2_bind[,"parameter"] %in% c("Nvalid", "Ncases", "R2", "R2nagel")),]
  }
  # drop significance
  jk2_bind <- jk2_bind[jk2_bind[["coefficient"]] %in% c("est", "se"), ]
  # reshape
  jk2_wide <- dcast ( jk2_bind, as.formula ( paste ( " ... ~ coefficient + ", tv ,sep="") ) )
  # calculate trend
  lev <- unique(jk2_bind[[tv]])
  jk2_wide[,"est_trend"] <- jk2_wide[,paste("est_",lev[2],sep="")] - jk2_wide[,paste("est_",lev[1],sep="")]
  # merge linking errors
  jk2_wide <- merge(jk2_wide, le, by = c("parameter", "depVar"), all = TRUE)
  # calculate trend SEs
  jk2_wide[,"se_trend"] <- sqrt(jk2_wide[, paste("se_",lev[2],sep="")]^2 + jk2_wide[, paste("se_",lev[1],sep="")]^2 + jk2_wide[, "le"]^2)
  # calculate p-Values
  jk2_wide[, "sig_trend"] <- 2*pnorm(abs(jk2_wide[, "est_trend"]/jk2_wide[, "se_trend"]), lower=FALSE)
  # Effect size for means (only for jk2.mean)
  es <- character(0)
  if (  fun == "mean" ) {
    jk2_wide[, "es_trend"]  <- NA
    # not for groupDiffs as no SD is provided by eatRep, split up data frame and rbind later
    jk2_wideS <- jk2_wide[!jk2_wide[, "comparison"] %in% c("crossDiff_of_groupDiff", "groupDiff") | is.na(jk2_wide[, "comparison"]), ]
    pooledSD <- sqrt(0.5 * (jk2_wideS[jk2_wideS[, "parameter"] == "sd", paste("est_",lev[1], sep="")]^2 +
                       jk2_wideS[jk2_wideS[, "parameter"] == "sd", paste("est_",lev[2], sep="")]^2))
    jk2_wideS[jk2_wideS[, "parameter"] == "mean", "es_trend"]  <- jk2_wideS[jk2_wideS[, "parameter"] == "mean", "est_trend"] / pooledSD
    es <- "es_trend" # add es to reshaping!
    jk2_wide <- unique(rbind(jk2_wide, jk2_wideS))
  }
  # reshape back to standard format
  jk2_add <- melt ( jk2_wide, measure.vars = c(paste("est_",lev[1], sep=""),
                                                         paste("est_",lev[2], sep=""),
                                                         paste("se_",lev[1], sep=""),
                                                         paste("se_",lev[2], sep=""),
                                                         "est_trend", "se_trend", "sig_trend", es), na.rm = FALSE)
  jk2_add <- jk2_add[!is.na(jk2_add[, "value"]), ] # drop missing rows (e.g. for SD, as no es can be calc)
  # split up variable column into coef and trend variable
  jk2_add <- separate(jk2_add, col = "variable", into = c("coefficient", tv), sep = "_")
  # rbind old results (drops le)
  jk2_out <- unique(rbind(jk2_add[, names(jk2_all)], jk2_all))
  return(jk2_out)
}
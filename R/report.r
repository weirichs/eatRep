report <- function ( repFunOut, trendDiffs = FALSE, add=list(), exclude = c("Ncases", "NcasesValid", "var", "sampleSize"), printGlm = FALSE, round = TRUE, digits = 3, printDeviance = FALSE) {
          if(is.null(repFunOut)) {return(NULL)}
          if ( length(grep("glm", as.character(repFunOut[["resT"]][[1]][1,"modus"]))) ==1 ) {
               if ( printGlm == TRUE ) { dG(repFunOut, digits = digits, printDeviance = printDeviance, add = add ) }
          }
          jk2      <- repFunOut[["resT"]]
          tv       <- repFunOut[["allNam"]][["trend"]]
          cols     <- c("group", "depVar", "modus", "parameter")
          grpv     <- setdiff(setdiff(colnames(jk2[[1]]), cols), c("comparison", "coefficient", "value", tv))
          grp_by   <- repFunOut[["allNam"]][["group.differences.by"]]
          cl_diffs <- repFunOut[["allNam"]][["cross.differences"]]
          funs     <- c("mean", "table", "quantile", "glm")
          fun      <- funs [ which( unlist(lapply(funs, FUN = function ( f ) { length(grep(f, jk2[[1]][1,"modus"]))})) > 0) ]
          
          if ( is.list(cl_diffs) ) {
               jk2 <- lapply(jk2, FUN = function (df) {computeCrossLevel (df, cols=cols, grpv = grpv, fun = fun, cl_diffs = cl_diffs, comp_type = "crossDiff")})
               jk2 <- lapply(jk2, FUN = function (df) {                         
                      spl <- strsplit(df[,"group"], ".vs.")
                      ind <- which(sapply(spl, length)==4)
                      if ( length(ind)>0) {
                           df[ind,"comparison"] <- "crossDiff_of_groupDiff"     ### das mittlere ".vs." gross schreiben
                           df[ind,"group"]      <- unlist(lapply(spl[ind], FUN = function ( z ) {paste(z[1], ".vs.", z[2], ".VS.", z[3], ".vs.", z[4], sep="")}))
                      }
                      return(df)})
          }
          if(!is.null(repFunOut[["SE_correction"]]) && !is.null(repFunOut[["SE_correction"]][[1]])) {
            if(length(which(jk2[[1]][["comparison"]] == "crossDiff_of_groupDiff")) > 0 ) {
              warning("Standard error correction for 'crossDiff_of_groupDiff' is currently not supported.")
            }
            mult_hierarchy <- any(unlist(lapply(repFunOut$allNam$cross.differences, function(x) x[2] - x[1] != 1)))
            if(mult_hierarchy) warning("Standard error correction for crossDifferences across multiple hierarchy levels is currently not supported.")

            jk2 <- lapply(jk2, function(jk2_single) {
              seCorrect(SE_correction = repFunOut[["SE_correction"]], jk2 = jk2_single, grpv = grpv)
            })
          }
          if ( !is.null(tv) ) {
               jk2 <- computeTrend(jk2 = jk2, le = repFunOut[["le"]], tv = tv, fun = fun)
          } else {
               jk2 <- jk2[[1]]
          }
          if ( !is.null(tv) && trendDiffs ) {
               jk2 <- computeTrendDiffs(jk2 = jk2, grpv = grpv, tv = tv, grp_by = grp_by, fun = fun, cl_diffs = cl_diffs)
          }
          if ( length(add)>0) {
               if(!all(nchar(names(add))>0)) { stop("'add' must be named.")}    
               if(length(names(add)) != length(unique(names(add)))) { stop("Duplicated names of 'add' are not allowed.")}
               if(!all(sapply(add, length) == 1)) {stop("All elements of 'add' must be of length 1.")}
               if(!all(sapply(add, class) == "character")) {stop("All elements of 'add' must be of class 'character'.")}
               dopp<- names(add) %in% colnames(jk2)
               ind <- which(dopp==TRUE)
               if ( length( ind ) > 0 ) {stop(paste0("Following names of 'add' are not allowed: '",paste(names(add)[ind], collapse = "', '"), "'."))}
               for ( u in names(add)) {jk2[,u] <- add[[u]]}
          }
          spltVar  <- c("coefficient", tv)                                      
          if ( length(exclude)>0) {
               weg <- which(jk2[,"parameter"] %in% exclude)
               if ( length(weg)>0) {
                    jk2 <- jk2[-weg,]
               }
          }                                                                     
          frml     <- as.formula(paste0("... ~ ", paste(spltVar,collapse=" + ") ) )
          jk2wide  <- dcast(data = jk2, formula = frml, value.var = "value")
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
                      y  <- z[1,]                                               
                      y[["coefficient"]] <- "p"                                 
                      y[["value"]]       <- 2*pnorm(abs(z[which(z[,"coefficient"] == "est"),"value"] / z[which(z[,"coefficient"] == "se"),"value"]), lower.tail=FALSE)
                      x  <- rbind ( x, y)                                       
                 }
                 return(x)}))                                                   
          return(dat)}


seCorrect <- function( SE_correction, jk2, grpv ) {
  UseMethod("seCorrect")
}

seCorrect.old <- function( SE_correction, jk2, grpv ) {
  jk2
}

seCorrect.wec_se_correction <- function( SE_correction, jk2, grpv ) {
  sep_jk2 <- separate_jk2(jk2 = jk2)
  cross_diff <- sep_jk2[["cross_diff"]]

  for(i in seq_along(SE_correction)) {
    output <- SE_correction[[i]][["resT"]][[sep_jk2[["year"]]]]

    rows <- length(SE_correction[[i]][["vgl"]][["groups.divided.by"]])
    single_grpv <- SE_correction[[i]][["focGrp"]]

    SEs <- output[!output$parameter %in% c("(Intercept)", "Nvalid", "R2"), c("parameter", "value", "coefficient")]
    SEs[, "parameter"] <- gsub(single_grpv, "", SEs[, "parameter"])
    SEs <- as.data.frame(pivot_wider(SEs, names_from = "coefficient", values_from = "value"))

    for(param in SEs[["parameter"]]) {
      esc_param <- escapeRegex(param)
      if(identical(SE_correction[[i]]$refGrp, "all")) { 
        grp_regexp <- paste0("^", esc_param, "\\.vs")
        compare_point_estimates(old_est = cross_diff[cross_diff$parameter == "mean" & grepl(grp_regexp, cross_diff$group) & cross_diff$coefficient == "est", "value"],
                                new_est = SEs[SEs[, "parameter"] == param, "est"],
                                param = param)
        cross_diff[cross_diff$parameter == "mean" & grepl(grp_regexp, cross_diff$group) & cross_diff$coefficient == "se",
                   "value"] <- SEs[SEs[, "parameter"] == param, "se"]
        cross_diff[cross_diff$parameter == "mean" & grepl(grp_regexp, cross_diff$group) & cross_diff$coefficient == "p",
                   "value"] <- SEs[SEs[, "parameter"] == param, "p"]
      } else { 


        param_selector <- paste0("^", esc_param, "\\.|", "^", esc_param, "_|", "_", esc_param, "\\.|", "_", esc_param, "_")
        col_names <- SE_correction[[i]]$refGrp[, "groupName"]
        col_levels <- SE_correction[[i]]$refGrp[, "groupValue"]
        filt_var <- recursive_filter(df = cross_diff, vars = col_names, var_levels = col_levels)
        compare_point_estimates(old_est = cross_diff[which(filt_var & grepl(param_selector, cross_diff$group) & cross_diff$coefficient == "est"), "value"],
                                new_est = SEs[SEs[, "parameter"] == param, "est"],
                                param = param)
        cross_diff[which(filt_var & grepl(param_selector, cross_diff$group) & cross_diff$coefficient == "se"),
                   "value"] <- SEs[SEs[, "parameter"] == param, "se"]
        cross_diff[which(filt_var & grepl(param_selector, cross_diff$group) & cross_diff$coefficient == "p"),
                   "value"] <- SEs[SEs[, "parameter"] == param, "p"]
      }
    }
  }
  cross_diff$modus <- unique(SE_correction[[1]][[1]][[1]]$modus)
  
  rbind(sep_jk2[["no_cross_diff"]], cross_diff)
}

separate_jk2 <- function(jk2) {
  year <- as.character(unique(jk2[["year"]]))
  if(length(year) == 0) year <- 1
  no_cross_diff <- jk2[is.na(jk2$comparison) | jk2$comparison != "crossDiff" | jk2$parameter != "mean", ]
  cross_diff <- jk2[which(jk2$comparison == "crossDiff" & jk2$parameter == "mean"), ]
  stopifnot(identical(nrow(no_cross_diff) + nrow(cross_diff), nrow(jk2)))
  stopifnot(nrow(cross_diff) > 0)
  list(cross_diff = cross_diff, no_cross_diff = no_cross_diff, year = year)
}

recursive_filter <- function(df, vars, var_levels) {
  stopifnot(length(vars) == length(var_levels))
  filt_var <- rep(TRUE, nrow(df))
  for(i in seq_along(vars)) {
    filt_var <- filt_var & df[[vars[i]]] == var_levels[[i]]
  }
  filt_var
}

seCorrect.rep_se_correction <- function( SE_correction, jk2, grpv ) {
  sep_jk2 <- separate_jk2(jk2 = jk2)
  cross_diff <- sep_jk2[["cross_diff"]]
  no_cross_diff <- sep_jk2[["no_cross_diff"]]
  for(i in seq_along(SE_correction)) {
    output <- SE_correction[[i]][["resT"]][[sep_jk2[["year"]]]]
    rows <- length(SE_correction[[i]][["vgl"]][["groups.divided.by"]])
    single_grpv <- SE_correction[[i]][["focGrp"]]
    for(param in output[["group"]]) {
      compare_point_estimates(old_est = cross_diff[cross_diff$group == param & cross_diff$coefficient == "est", "value"],
                              new_est = output[output$group == param & output$coefficient == "est", "value"],
                              param = param)
      cross_diff[cross_diff$group == param & cross_diff$coefficient == "se", "value"] <- output[output$group == param & output$coefficient == "se", "value"]
      cross_diff[cross_diff$group == param & cross_diff$coefficient == "p", "value"]  <- output[output$group == param & output$coefficient == "p", "value"]
    }
  }
  rbind(no_cross_diff, cross_diff)
}

compare_point_estimates <- function(old_est, new_est, param) {
  if(abs(abs(old_est) - abs(new_est)) >= 1e-10) {
    warning("Difference in point estimate of cross level difference for comparison ", param, ": ", round(abs(old_est) - abs(new_est), digits = 3))
  }
  return()
}

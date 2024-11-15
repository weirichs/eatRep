report2 <- function ( repFunOut, add=list(), exclude = c("NcasesValid", "var"), printGlm = FALSE,
                     round = TRUE, digits = 3, printDeviance = FALSE, printSE_correction = FALSE) {
          if(is.null(repFunOut)) {return(NULL)}
          allN     <- repFunOut[["allNam"]]
          repFunOut<- buildList(repFunOut)
          out      <- do.call("rbind", lapply(names(repFunOut), FUN = function (rfo) {
                      if ( length(grep("glm", as.character(repFunOut[[rfo]][["resT"]][[1]][1,"modus"]))) ==1 ) {
                           if ( printGlm == TRUE ) { dG(repFunOut[[rfo]], digits = digits, printDeviance = printDeviance, add = add ) }
                      }
                      jk2      <- repFunOut[[rfo]][["resT"]]
                      tv       <- repFunOut[[rfo]][["allNam"]][["trend"]]
                      cols     <- c("group", "depVar", "modus", "parameter")
                      grpv     <- setdiff(setdiff(colnames(jk2[[1]]), cols), c("comparison", "coefficient", "value", tv))
                      grp_by   <- repFunOut[[rfo]][["allNam"]][["group.differences.by"]]
                      if(is.logical(repFunOut[[rfo]][["allNam"]][["cross.differences"]]) &&  isTRUE(repFunOut[[rfo]][["allNam"]][["cross.differences"]])) {
                        cl_diffs <- combinat::combn(0:length(repFunOut[[rfo]][["allNam"]][["group"]]),2, simplify=FALSE)
                      }  else  {
                        cl_diffs <- repFunOut[[rfo]][["allNam"]][["cross.differences"]]
                      }
                      funs     <- c("mean", "table", "quantile", "glm", "lmer")
                      fun      <- funs [ which( unlist(lapply(funs, FUN = function ( f ) { length(grep(f, jk2[[1]][1,"modus"]))})) > 0) ]
                      if ( fun == "table") {jk2 <- lapply(jk2, reduceDoubleN)}
                      if ( is.list(cl_diffs) ) {
                           jk2 <- lapply(jk2, FUN = function (df) {
                                  ret <- plyr::rbind.fill ( df, computeCrossLevel (df, cols=cols, grpv = grpv, fun = fun, cl_diffs = cl_diffs, allNams = allN))
                                  ret[,"row"] <- 1:nrow(ret)                    
                                  return(ret) })
                      }
                      if(!is.null(repFunOut[[rfo]][["SE_correction"]]) && !is.null(repFunOut[[rfo]][["SE_correction"]][[1]])) {
                          if(length(which(jk2[[1]][["comparison"]] == "crossDiff_of_groupDiff")) > 0 ) {
                              warning("Standard error correction for 'crossDiff_of_groupDiff' is currently not supported.")
                          }
                          mult_hierarchy <- any(unlist(lapply(cl_diffs, function(x) x[2] - x[1] != 1)))
                          if(mult_hierarchy) {warning("Standard error correction for crossDifferences across multiple hierarchy levels is currently not supported.")}
                          jk2 <- lapply(jk2, function(jk2_single) { seCorrect(SE_correction = repFunOut[[rfo]][["SE_correction"]], jk2 = jk2_single, grpv = grpv, allNam=allN, printSE_correction=printSE_correction) })
                      }
                      if ( !is.null(tv) ) {
                           jk2 <- computeTrend(jk2 = jk2, repFunOut = repFunOut[[rfo]], tv = tv, fun = fun, allNam=allN)
                      } else {
                           jk2 <- jk2[[1]]
                      }
                      if ( length(add)>0) {                                     
                           checkmate::assert_list(add, types="character", unique=TRUE, min.len = 1, names = "unique")
                           if(!all(sapply(add, length) == 1)) {stop("All elements of 'add' must be of length 1.")}
                           dopp<- names(add) %in% colnames(jk2)
                           ind <- which(dopp==TRUE)
                           if ( length( ind ) > 0 ) {stop(paste0("Following names of 'add' are not allowed: '",paste(names(add)[ind], collapse = "', '"), "'."))}
                           for ( u in names(add)) {jk2[,u] <- add[[u]]}
                      }
                      jk2wide  <- reshape2::dcast(data = jk2[,-eatTools::whereAre(c("group","row", "hierarchy.level"), colnames(jk2), verbose=FALSE)], formula = ... ~ coefficient, value.var = "value")
                      if ( fun == "glm") {
                           jk2wide[,"parameter"] <- car::recode(jk2wide[,"parameter"], "'Nvalid'='zzzzNvalid'; 'R2'='zzzzR2'; 'R2nagel'='zzzzR2nagel'")
                           jk2wide <- data.frame(jk2wide[sort(jk2wide[,"parameter"],decreasing=FALSE,index.return=TRUE)$ix,])
                           jk2wide[,"parameter"] <- car::recode(jk2wide[,"parameter"], "'zzzzNvalid'='Nvalid'; 'zzzzR2'='R2'; 'zzzzR2nagel'='R2nagel'")
                      }
                      if ( isTRUE(round)) {
                           jk2wide <- eatTools::roundDF(jk2wide, digits = digits)
                      }
                      return(jk2wide)}))
          cols     <- grep("^id$|^unit_1$|^unit_2$", colnames(out), value=TRUE)
          alt      <- unique(unlist(out[,cols]))
          altneu   <- data.frame (alt=alt, neu = as.character(as.numeric(as.factor(alt))), stringsAsFactors = FALSE)
          out[,"id"] <- paste(car::recode(out[,"comparison"], "'none'='group'; else='comp'"),eatTools::recodeLookup(out[,"id"], altneu), sep="_")
          if ( all(c("unit_1", "unit_2") %in% colnames(out))) {
              for ( i in 1:2) {
                   out[,paste("unit", i, sep="_")] <- car::recode(paste(car::recode(out[,"comparison"], "'crossDiff'='group'; 'groupDiff'='group'; 'trend'='group'; 'none'='NA'; else='comp'"),eatTools::recodeLookup(out[,paste("unit", i, sep="_")], altneu), sep="_"), "'NA_NA'=NA")
              }
          }
          plain    <- data.frame ( label1 = createLabel1(out, allNam=allN), label2 = createLabel2(out, allNam=allN), out, stringsAsFactors=FALSE) |> dplyr::select(-dplyr::any_of(c("type")))
          if ( length(exclude)>0) {
               weg <- which(plain[,"parameter"] %in% exclude)
               if ( length(weg)>0) {plain <- plain[-weg,]}
          }
          compar   <- subset(plain, comparison!="none") |> dplyr::select(dplyr::any_of(c("id", "unit_1", "unit_2", "comparison"))) |> unique()
          groups   <- subset(plain, comparison=="none") |> dplyr::select(dplyr::any_of(c("id", allN[["group"]], allN[["trend"]],  names(add)))) |> unique()
          estim    <- plain |> dplyr::select(dplyr::any_of(c("id", "depVar", "parameter", "est", "se", "p", "es")))
          rownames(plain) <- NULL
          ret      <- list(plain = plain, comparisons = compar, group=groups, estimate = estim)
          class(ret) <- c("list", "report2")
          return(ret)}

computeCrossLevel <- function ( jk2, cols, grpv, fun, cl_diffs, allNams) {
       ori <- jk2
       if ( fun == "mean" ) {
            jk2  <- jk2[which(jk2[,"parameter"] %in% c("mean", "sd")),]
       }
       if ( fun %in% c("glm", "table") ) {
            jk2  <- jk2[which(!jk2[,"parameter"] %in% c("Nvalid", "Ncases", "R2", "R2nagel")),]
       }
       if(length(allNams[["group"]]) < 2) {jk2 <- jk2[which(jk2[,"comparison"] == "none"),]}
       stopifnot(length(unique(jk2[,"comparison"])) <= 2)
       ret <- do.call("rbind", by ( data = jk2, INDICES = jk2[,"comparison"], FUN = function ( d ) {
              if ( d[1,"comparison"] != "none") { if(inherits(try(cl_diffs  <- combinat::combn(unique(d[,"hierarchy.level"]),2, simplify=FALSE)  ),"try-error"))  {return(NULL)}}
              ret2 <- do.call("rbind", lapply ( cl_diffs, FUN = function ( comp_vec ) {
                      redDF_1 <- d[which(d[, "hierarchy.level"] %in% comp_vec),]
                      fac     <- by(data = redDF_1, INDICES = redDF_1[, "hierarchy.level"], FUN = function ( x ) { unique(x[, "group"])})
                      if(length(fac) == 1) {return(NULL)}
                      all_grp_list <- do.call("rbind", lapply(fac[[1]], function(high_lvl) {
                          hl_levels <- unlist(strsplit(high_lvl, ", |_+"))
                          matches   <- unlist(lapply(fac[[2]], FUN = function (single) {
                                       part <- unlist(strsplit(single, ", |_+"))
                                       return(sum(part %in% hl_levels))}))
                          low_lvl   <- fac[[2]][which(matches == max(matches))]
                          vgl       <- expand.grid(high_lvl, low_lvl)
                          grp       <- do.call("rbind", plyr::alply(as.matrix(vgl), .margins = 1, .fun = function(single_comp) {
                                       redDF_2 <- redDF_1[which(redDF_1[,"group"] %in% single_comp),]
                                       redDF_2 <- redDF_2[c(which(redDF_2[, "hierarchy.level"] == comp_vec[1]), which(redDF_2[, "hierarchy.level"] == comp_vec[2])), ]
                                       newRows <- compareParameters(df_allP = redDF_2, grpv = grpv, fun = fun, allNams=allNams)
                                       return(newRows)  }))
                          return(grp) }))
                      return(all_grp_list)}))
              return(ret2)}))
       return(ret)}

compareParameters <- function(df_allP, grpv, fun, allNams) {
  df_allP <- df_allP[with(df_allP, order(parameter, group)),]
  df_allP <- df_allP[which(df_allP[,"coefficient"] %in% c("est", "se")),]       
  if(identical(fun, "mean")) { df_sd <- df_allP[df_allP[, "parameter"] == "sd", ] }
  out <- do.call("rbind", by ( data = df_allP, INDICES = df_allP[,"parameter"], FUN = function ( df ) {
     df   <- cleanDF(df)                                                        
     df   <- df[order(df[,"hierarchy.level"], decreasing = FALSE), ]            
     meaD <- diff(df[which(df[,"coefficient"] == "est"),"value"])               
     if ( length(meaD) ==0 || is.na(meaD)) {return(NULL)}
     if ( !"se" %in% df[, "coefficient"]) {
         warning( "No standard error for parameter '",unique(df[,"parameter"]),"'. Cannot compute standard errors and p value for cross-level difference between '",df[1,"group"],"' and '",df[2,"group"],"'.\n")
         se <- pval <- NA
     }  else  {
         se  <- sqrt(sum(df[which(df[,"coefficient"] == "se"),"value"]^2))      
         pval<- 2*pnorm(abs(meaD/se), lower.tail=FALSE)
     }
     es  <- NA                                                                  
     if ( fun == "mean" && df[1,"parameter"] == "mean" && nrow(df_sd)>0 ) {
         sd_wide <- reshape2::dcast(df_sd[which(df_sd[,"coefficient"] == "est"),], group~parameter, value.var = "value")
         if ( nrow(sd_wide) > 1) { es  <- meaD / sqrt(0.5*sum(sd_wide[,"sd"]^2)) }
     }
     if(nrow(df) == 2 ) {                                                       
       ret <- rbind (df, df)                                                    
     }  else  {
       ret <- df
     }
     ret[,"comparison"]  <- car::recode(ret[,"comparison"], "'none'='crossDiff'; 'groupDiff'='crossDiff_of_groupDiff'")
     ret[,"coefficient"] <- c("est","se", "p", "es")
     ret[,"value"]       <- c(meaD, se, pval, es)
     const <- unlist(lapply(allNams[["group"]], FUN = function (v) {length(unique(ret[,v])) ==1}))
     vars  <- allNams[["group"]][which(!const)]
     foo   <- lapply(vars, FUN = function (v) {stopifnot(any(grepl("total", ret[,v])))})
     for (v in vars) {ret[,v] <- paste0(unique(setdiff(ret[,v], "total")),  " - total")}
     ret[,"unit_1"] <- sort(unique(ret[,"id"]))[1]
     ret[,"unit_2"] <- sort(unique(ret[,"id"]))[2]
     ret[,"id"]     <- paste(genTS(), paste(as.character(ret[,"row"]), collapse=""), sep="_")
     ret   <- eatTools::na_omit_selection(ret, "value")
     return(ret)}))
  return(out)}

cleanDF <- function(df){
   if(nrow(df) %in% c(2,4)) {return(df)}                                        


   stopifnot(nrow(df) %in% c(1,3))
   valid <- table(df[,"coefficient"])
   weg   <- which(valid != 2)
   if ( length(weg) >0) {
         warning("No '",paste(names(weg), collapse = "', '"), "' for parameter '",unique(df[,"parameter"]),"'. Skip computation of standard errors and p values for cross-level difference between '",df[1,"group"],"' and '",df[2,"group"],"'")
         df <- df[-eatTools::whereAre(names(weg), df[,"coefficient"], verbose=FALSE),]
   }
   return(df)}


addSig <- function ( dat , groupCols = NULL , allNam = NULL ) {
  checkmate::assert_data_frame(dat)
  checkmate::assert_character(groupCols, null.ok = TRUE)
  
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


seCorrect <- function( SE_correction, jk2, grpv, allNam, printSE_correction ) {
      if(inherits(SE_correction, "old")) {return(jk2)}
      stopifnot(nrow(jk2) == length(unique(jk2[,"row"])))
      stopifnot(length(which(is.na(jk2[,"row"]))) == 0)
      filt       <- intersect(which(jk2[,"comparison"] == "crossDiff"), which(jk2[,"parameter"] == "mean"))
      cross_diff <- jk2[filt,]
      noCrossDiff<- jk2[setdiff(1:nrow(jk2), filt),]
      if ( !is.null(allNam[["trend"]])) {year <- as.character(unique(jk2[,allNam[["trend"]]]))} else {year <- 1}
      ret        <- seCorrect.wec_se_correction (SE_correction, jk2, grpv, allNam, cross_diff, noCrossDiff , year, printSE_correction)
      return(ret)}

seCorrect.wec_se_correction <- function( SE_correction, jk2, grpv, allNam, cross_diff, noCrossDiff, year, printSE_correction ) {
      for(i in seq_along(SE_correction)) {
          single_grpv<- SE_correction[[i]][["focGrp"]]
          output     <- SE_correction[[i]][["resT"]][[year]]
          if(inherits(SE_correction, "rep_se_correction")) {output[,"parameter"] <- output[,single_grpv]}
          rows       <- length(SE_correction[[i]][["vgl"]][["groups.divided.by"]])
          SEs        <- output[!output$parameter %in% c("(Intercept)", "Nvalid", "R2"), c("parameter", "value", "coefficient")]
          if(inherits(SE_correction, "wec_se_correction")) {SEs[, "parameter"] <- gsub(paste0("^", single_grpv), "", SEs[, "parameter"])}
          SEs        <- eatTools::makeDataFrame(tidyr::pivot_wider(SEs, names_from = "coefficient", values_from = "value"), verbose=FALSE)
          for(param in SEs[["parameter"]]) {
              esc_param <- Hmisc::escapeRegex(param)
              if(identical(SE_correction[[i]][["refGrp"]], "all")) {            
                  olds <- lapply(c("est", "se", "p"), FUN = function (coeff) {  
                          old  <- cross_diff[cross_diff$parameter == "mean" & cross_diff[,single_grpv] == paste0(esc_param," - total") & cross_diff$coefficient == coeff, ]
                          cols <- setdiff(allNam[["group"]], single_grpv)
                          if ( length(cols)>0) {
                               for ( co in cols) {old[,co] <- car::recode(old[,co], "'total'='total'; else=NA")}
                               old  <- eatTools::na_omit_selection(old, cols)
                          }
                          return(old)})
                  names(olds) <- c("est", "se", "p")
                  compare_point_estimates(old_est = olds[["est"]], new_est = SEs[which(SEs[, "parameter"] == param), ],old_col ="value", new_col = "est", param = param)
                  if(printSE_correction){cat(paste0(strsplit(class(SE_correction)[1], "_")[[1]][1],": Replace old SE ", paste(round(cross_diff[match(olds[["se"]][,"row"], cross_diff[,"row"]), "value"], digits = 4), collapse=", "), " with ", paste(round(SEs[SEs[, "parameter"] == param, "se"], digits = 4), collapse=", "), "\n"))}
                  cross_diff[match(olds[["se"]][,"row"], cross_diff[,"row"]), "value"] <- SEs[SEs[, "parameter"] == param, "se"]
                  if(printSE_correction){cat(paste0(strsplit(class(SE_correction)[1], "_")[[1]][1],": Replace old p value ", paste(round(cross_diff[match(olds[["p"]][,"row"], cross_diff[,"row"]), "value"], digits = 4), collapse=", "), " with ", paste(round(SEs[SEs[, "parameter"] == param, "p"], digits = 4), collapse=", "), "\n"))}
                  cross_diff[match(olds[["p"]][,"row"], cross_diff[,"row"]), "value"] <- SEs[SEs[, "parameter"] == param, "p"]
              } else {                                                          
                  old  <- cross_diff                                            
                  ref  <- SE_correction[[i]][["refGrp"]]
                  for ( reihe in 1:nrow(ref)) { old <- old[which(old[,ref[reihe,"groupName"]] == ref[reihe,"groupValue"]),]}
                  old  <- old[which(old[,SE_correction[[i]][["focGrp"]]] == paste0(esc_param, " - total")),]
                  compare_point_estimates(old_est = old[which(old[,"coefficient"] == "est"),], new_est = SEs[which(SEs[, "parameter"] == param), ],old_col ="value", new_col = "est", param = param)
                  if(printSE_correction){cat(paste0(strsplit(class(SE_correction)[1], "_")[[1]][1],": Replace old SE ", paste(round(cross_diff[match(old[which(old[,"coefficient"] == "se"),"row"], cross_diff[,"row"]), "value"], digits = 4), collapse=", "), " with ", paste(round(SEs[SEs[, "parameter"] == param, "se"], digits = 4), collapse=", "), "\n"))}
                  cross_diff[match(old[which(old[,"coefficient"] == "se"),"row"], cross_diff[,"row"]), "value"] <- SEs[SEs[, "parameter"] == param, "se"]
                  if(printSE_correction){cat(paste0(strsplit(class(SE_correction)[1], "_")[[1]][1],": Replace old p ", paste(round(cross_diff[match(old[which(old[,"coefficient"] == "p"),"row"], cross_diff[,"row"]), "value"], digits = 4), collapse=", "), " with ", paste(round(SEs[SEs[, "parameter"] == param, "p"], digits = 4), collapse=", "), "\n"))}
                  cross_diff[match(old[which(old[,"coefficient"] == "p"),"row"], cross_diff[,"row"]), "value"] <- SEs[SEs[, "parameter"] == param, "p"]
              }
          }
      }
return(rbind(noCrossDiff, cross_diff))}


compare_point_estimates <- function(old_est, new_est, param, old_col, new_col) {
  if(nrow(old_est) != nrow(new_est)) {
     rownames(old_est) <- rownames(old_est) <- NULL; cat(paste0("\n\nWarning: Length of old and new values differ. \n\n   Old: \n", eatTools::print_and_capture(old_est, 5), "\n\n   New: \n", eatTools::print_and_capture(new_est, 5)))
  }
  if(abs(abs(old_est[,old_col]) - abs(new_est[,new_col]))[1] >= 1e-6) {
     warning("Difference in point estimate of cross level difference for comparison ", param, ": ", round(abs(old_est[,old_col]) - abs(new_est[,new_col]), digits = 3))
  }
  return()}

reduceDoubleN <- function(jk2){
          cases1<- which(!duplicated(jk2[ ,-eatTools::whereAre(c("unit_1", "unit_2", "id"),colnames(jk2), verbose=FALSE)]))
          cases <- unique(jk2[intersect(which(jk2[,"parameter"] %in% c("Ncases", "NcasesValid")),cases1),])
          jk2   <- rbind(jk2[which(jk2[,"parameter"]  %nin% c("Ncases", "NcasesValid")),], cases)
          return(jk2)}

buildList <- function(repFunOut){
      dvs <- unique(repFunOut[["resT"]][[1]][,"depVar"])
      out <- lapply(dvs, FUN = function (dv) {
             for ( i in 1:length(repFunOut[["resT"]])) {repFunOut[["resT"]][[i]] <- repFunOut[["resT"]][[i]][which(repFunOut[["resT"]][[i]][,"depVar"] == dv),]}
             repFunOut[["allNam"]][["dependent"]] <- dv
             return(repFunOut)})
      names(out) <- dvs
      return(out)}

computeTrend <- function(jk2, tv, repFunOut, fun, allNam) {
        jk2_bind<- do.call("rbind", jk2)                                        
        if(identical(fun, "mean")) {                                            
            jk2_bind[,"splitVar"] <- jk2_bind[["parameter"]] %in% c("mean", "sd")
        }
        if ( fun %in% c("glm", "table") ) {
            jk2_bind[,"splitVar"] <- !jk2_bind[,"parameter"] %in% c("Nvalid", "Ncases", "R2")
        }
        if ( identical(fun, "lmer") ) {
            jk2_bind[,"splitVar"] <- !jk2_bind[,"parameter"] %in% c("Nvalid", "R2_Lev2","R2_Lev1","R2_Total","ICC_Uncond","ICC_UncondWB", "ICC_Cond")
        }                                                                       
        if ( identical(fun, "quantile") ) {
            jk2_bind[,"splitVar"] <- TRUE
        }
        jk2_bind[setdiff(1:nrow(jk2_bind), which(jk2_bind[,"coefficient"] %in% c("est", "se"))),"splitVar"] <- FALSE
        lev     <- unique(jk2_bind[,tv])                                        
        le      <- check2(repFunOut=repFunOut, jk2=jk2_bind, fun=fun, lev=lev)
        vgl <- combinat::combn(names(jk2),2, simplify=FALSE)                    
        adds<- unique(do.call("rbind", lapply( 1:length(vgl), FUN = function ( comp) {
               jk2_binS<- check1(jk2=jk2[vgl[[comp]]], jk2_bind=jk2_bind[intersect(which(jk2_bind[,"splitVar"] == TRUE), which(jk2_bind[,tv] %in% vgl[[comp]])),], tv=tv, allNam = allNam)
               wide    <- eatTools::makeDataFrame(tidyr::pivot_wider(jk2_binS[,-eatTools::whereAre(c("row","unit_1", "unit_2", "variable"), colnames(jk2_binS), verbose=FALSE)], names_from = c("coefficient",allNam[["trend"]]), values_from = c("value", "id")), verbose=FALSE)
               wide[,"est"]        <- wide[,paste0("value_est_",sort(vgl[[comp]])[2])] - wide[,paste0("value_est_",sort(vgl[[comp]])[1])]
               wide[,"comparison"] <- eatTools::removePattern(paste0("trend_",wide[,"comparison"]),"_none")
               ind     <- intersect(which(le[,"trendLevel1"] %in% vgl[[comp]]), which(le[,"trendLevel2"] %in% vgl[[comp]]))
               stopifnot(length(ind)>0)
               le_S    <- le[ind,]
               stopifnot(all(as.vector(unlist(by(le_S, INDICES = le_S[,c("parameter", "depVar")], nrow))) == 1))
               wide    <- merge(wide, le_S[,-na.omit(match(c("trendLevel1", "trendLevel2", "domain"), colnames(le_S)))], by = c("parameter", "depVar"), all = TRUE)
               miss    <- which(is.na(wide[,"le"]))
               if ( length(miss)>0){
                    warning(paste0("Found ",length(miss)," missing linking errors for dependent variable '",unique(wide[,"depVar"]),"' and parameter(s) '",paste(unique(wide[which(is.na(wide[,"le"])),"parameter"]), collapse="', '"),"'. Assume linking error of 0 for these cases."))
                    wide[which(is.na(wide[,"le"])),"le"] <- 0
               }
               wide[,"se"] <- sqrt(wide[, paste("value_se_",vgl[[comp]][1],sep="")]^2 + wide[, paste("value_se_",vgl[[comp]][2],sep="")]^2 + wide[, "le"]^2)
               wide[,"p"]  <- 2*pnorm(abs(wide[,"est"]/wide[, "se"]), lower.tail=FALSE)
               existSD <- "sd" %in% jk2_binS[,"parameter"]                      
               if(fun == "mean" && !existSD && comp == 1) {message("Cannot find standard deviations in output. Skip computation of effect sizes.")}
               if (  fun == "mean" && existSD) {                                
                   wide2  <- data.frame (parameter = "mean", tidyr::pivot_wider(jk2_binS[which(jk2_binS[,"comparison"] == "none"),-eatTools::whereAre(c("row","unit_1", "unit_2", "variable"), colnames(jk2_binS), verbose=FALSE)], names_from = c("parameter", "coefficient",allNam[["trend"]]), values_from = c("value", "id")), stringsAsFactors = FALSE)
                   altneu <- data.frame (  alt = paste0("id_mean_est_", vgl[[comp]][1]), neu = paste0("id_est_", vgl[[comp]][1]), stringsAsFactors = FALSE)
                   colnames(wide2) <- eatTools::recodeLookup(colnames(wide2), altneu)
                   stopifnot(length(unique(wide2[, altneu[["neu"]]])) == nrow(wide2) )
                   wide2[,"pooledSD"] <- sqrt(0.5 * (wide2[,paste0("value_sd_est_",vgl[[comp]][1])]^2 + wide2[,paste0("value_sd_est_",vgl[[comp]][2])]^2))
                   wide2[,"es"]       <- (wide2[,paste0("value_mean_est_",sort(vgl[[comp]])[2])] - wide2[,paste0("value_mean_est_",sort(vgl[[comp]])[1])]) / wide2[,"pooledSD"]
                   wide   <- eatTools::mergeAttr(wide,wide2[,intersect(c(allNam[["group"]], altneu[["neu"]],"parameter", "es"), colnames(wide2))], all=TRUE, setAttr=FALSE, xName = "original", yName = "effectSize", verbose=c("unique", "common"))
               }                                                                
               mvs     <- intersect(c(allNam[["group"]], "est", "se", "p", "es"), colnames(wide))
               jk2Long <- reshape2::melt(data.frame ( wide[,c(paste0("id_est_",vgl[[comp]][1]) , paste0("id_est_",vgl[[comp]][2]), "depVar","group", "modus", "parameter", "comparison",mvs)], matrix (paste0(vgl[[comp]][2], " - ", vgl[[comp]][1]), ncol = 1, nrow = nrow(wide), dimnames = list(NULL, allNam[["trend"]])),stringsAsFactors = FALSE), measure.vars = setdiff(mvs,allNam[["group"]]), na.rm=TRUE, variable.name = "coefficient")
               colnames(jk2Long) <- eatTools::recodeLookup(colnames(jk2Long), data.frame ( alt = c(paste0("id_est_",vgl[[comp]][1]) , paste0("id_est_",vgl[[comp]][2])), neu = c("unit_1", "unit_2"), stringsAsFactors = FALSE))
               return(jk2Long)})))
        adds[,"id"] <- paste(adds[,"unit_1"], adds[,"unit_2"],sep="_")
        adds<- eatTools::rbind_common(adds, jk2_bind)
        return(adds) }                                                          

check1 <- function(jk2, jk2_bind, tv, allNam) {
       stopifnot(length(jk2) == 2)                                              
       mv  <- intersect(c("type", allNam[["group"]], "depVar", "modus", "comparison", "parameter", "coefficient"), colnames(jk2[[1]]))
       dupl<- lapply(jk2, FUN = function (y) {which(!duplicated(y[,c(mv,tv)]))})
       mrge<- merge(jk2[[1]][dupl[[1]],c(mv,tv, "row")], jk2[[2]][dupl[[2]],c(mv,tv, "row")], by = mv, all=FALSE, suffixes = names(jk2))
       stopifnot(nrow(mrge) <= max(sapply(jk2, nrow)))
       oriN<- which(sapply(names(jk2), FUN = function (x) {nrow(jk2[[x]][dupl[[x]],]) > nrow(mrge)}))
       if(length(oriN)>0) {
           for ( nam in names(oriN)) {
                weg <- setdiff(jk2[[nam]][,"row"], mrge[,paste0("row",nam)])
                message(paste0("\n   Following ",length(weg)," units in trend group '",nam, "' without counterpart in trend group '",setdiff(names(jk2), nam),"'.\n\n",eatTools::print_and_capture(jk2[[nam]][which(jk2[[nam]][,"row"] %in% weg),c("group", "depVar", "comparison", "parameter", "coefficient", "value")], spaces = 8), "\n\n   No '",paste(names(jk2), collapse = ".vs."),"' trends will be computed."))
           }
           mrgL<- reshape2::melt(mrge, id.vars = intersect(c(allNam[["group"]], "comparison", "parameter", "coefficient"),colnames(mrge)), measure.vars = grep("^year", colnames(mrge), value=TRUE), value.name = "year")
           jk2B<- eatTools::mergeAttr(jk2_bind, mrgL, all=FALSE, setAttr=FALSE)
       }  else  {
           jk2B<- jk2_bind
       }
       return(jk2B)}

check2 <- function(repFunOut, jk2, fun, lev){
       wdf <- attr(repFunOut[["le"]], "linkingErrorFrame")                      
       if(is.null(wdf)) {
           return(repFunOut[["le"]])
       }  else  {
           allV  <- list(trendLevel1 = "trendLevel1", trendLevel2 = "trendLevel2", parameter = "parameter", linkingError="linkingError", depVar = "depVar")
           allN  <- lapply(allV, FUN=function(ii) {eatTools::existsBackgroundVariables(dat = repFunOut[["le"]], variable=ii, warnIfMissing = TRUE)})
           le    <- repFunOut[["le"]]
           dv    <- unique(jk2[,"depVar"])
           stopifnot(length(dv)==1)
           if(!dv %in% le[,"depVar"]) {                                         
               le[,"depVar"] <- paste0(le[,"depVar"], le[,"parameter"])
               if(!dv %in% le[,"depVar"]) {stop(paste0("Cannot found dependent variable '",dv,"' in 'depVar' column of linking error data.frame 'linkErr'."))}
           }
           le  <- le[which(le[,"depVar"] == dv),]
           le1   <- le[,c("parameter", "trendLevel1", "trendLevel2")]
           if(nrow(le1) != nrow(unique(le1))) {stop("Linking error data.frame 'linkErr' is not unique in 'parameter', 'trendLevel1', 'trendLevel2'. Probable cause: 'linkErr' contains linking errors of multiple competence domains?")}
           add <- setdiff(unique(jk2[,"parameter"]), unique(le[,"parameter"]))
           if ( length(add)>0) { warning(paste0("No linking errors for parameters '",paste(add, collapse="', '"),"'. Linking errors for these parameters will be defaulted to 0."))}
           years <- eatTools::asNumericIfPossible(sort(unique (jk2[,repFunOut[["allNam"]][["trend"]]])), force.string=FALSE)
           combs <- combinat::combn(x=years, m=2, simplify=FALSE)
           exist <- suppressWarnings(unique(plyr::alply(eatTools::asNumericIfPossible(le,force.string=FALSE), .margins = 1, .fun = function (zeile) {sort(c(zeile[["trendLevel1"]], zeile[["trendLevel2"]]))})))
           drin  <- combs %in% exist
           if(!all(drin)) {stop(paste0("Data contains combination of years '",paste(combs[!drin], collapse="', '"), "' which are not included in linking error data.frame 'linkErr'."))}
           colnames(le) <- car::recode(colnames(le), "'linkingError'='le'")
           return(le)
       }}

createLabel1 <- function(plain, allNam) {
           label1 <- apply(X=plain, MARGIN = 1, FUN = function (zeile){
                     if(zeile[["parameter"]] == "chiSquareTest") {return("chiSquareTest")}
                     if(zeile[["comparison"]] == "none") {
                        if( length(allNam[["group"]])>0 && !identical (allNam[["group"]], "wholeGroup") ) {
                           no <- paste(zeile[allNam[["group"]]], collapse="_")
                        } else {
                           no <- "total"
                        }
                     } else {
                        no <- ""
                     }
                     if(!is.null(allNam[["trend"]]) && grepl(" - ", zeile[[allNam[["trend"]]]])) {
                        tr <- paste0("trend (",zeile[[allNam[["trend"]]]],") for ")
                        if(zeile[["comparison"]] == "trend") {
                           if(length(allNam[["group"]])>0 && !identical (allNam[["group"]], "wholeGroup")  ) {
                              tr <- paste0(tr, paste(zeile[allNam[["group"]]], collapse="_"))
                           } else {
                              tr <- paste0(tr, "total")
                           }
                        }
                     } else {
                        tr <- ""
                     }
                     if(grepl("crossdiff", zeile[["comparison"]],ignore.case=TRUE)) {
                        cdv<- setdiff(allNam[["group"]], allNam[["group.differences.by"]])
                        if ( length(cdv)>0) {
                            cd <- lapply(cdv, FUN = function (v) {
                                    spl <- unlist(strsplit(zeile[[v]], " - "))
                                    if(length(spl)==1) {spl <- c(spl,spl)}
                                    return(spl)})
                            cd1<- paste(unlist(lapply(cd, FUN = function (x) {x[1]})), collapse="_")
                            cd2<- paste(unlist(lapply(cd, FUN = function (x) {x[2]})), collapse="_")
                            cd <- paste0("crossDiff (",paste(cd1, cd2, sep=" - "),") ")
                        }  else  {
                            cd <- paste0("crossDiff (",zeile[[allNam[["group.differences.by"]]]], ")")
                        }
                     } else {
                        cd <- ""
                        tr <- eatTools::crop(tr, "for ")
                     }
                     if(grepl("groupdiff", zeile[["comparison"]],ignore.case=TRUE)) {
                        sig<- unlist(zeile[allNam[["group"]]])
                        col<- setdiff(grep(" - ", sig), grep("total", sig))
                        stopifnot(length(col)==1)
                        gd <- paste0("of groupDiff (",sig[[col]],")")           
                        ogc<- setdiff(allNam[["group"]], allNam[["group.differences.by"]])
                        ogc<- unlist(zeile[ogc])
                        weg<- c(grep(" - ", ogc), grep("total", ogc))
                        if(length(weg)>0) {ogc <- ogc[-weg]}
                        if(length(ogc)>0) {
                           gd <- paste0(gd, " in ", paste(names(ogc), ogc, sep="=", collapse=", "))
                        }
                        if(!is.null(allNam[["trend"]]) && length(grep(" - ", zeile[[allNam[["trend"]]]])) ==0) {
                           gd <- paste0(gd, " for ", allNam[["trend"]], " ", zeile[[allNam[["trend"]]]])
                        }
                     } else {
                        gd <- ""
                     }
                     if (tr=="" && cd =="") {gd <- eatTools::crop(gd, "of ") }
                     tot <- paste0(no, tr, cd, gd)
                     if (!is.null(allNam[["trend"]]) && length(grep(" - ", zeile[[allNam[["trend"]]]])) ==0) {
                         tot <- paste0(tot, " for ",allNam[["trend"]]," ",zeile[[allNam[["trend"]]]])
                     }
                     return(tot)})
           return(label1)}

createLabel2 <- function(plain, allNam) {
           gv     <- allNam[["group.differences.by"]]
           tv     <- allNam[["trend"]]
           label2 <- apply(X=plain, MARGIN = 1, FUN = function (zeile){
                     if(zeile[["parameter"]] == "chiSquareTest") {return("chiSquareTest")}
                     if(zeile[["comparison"]] == "none") {                      
                        if(!is.null(tv)) {
                           pre <- paste0(tv,"=",zeile[[tv]],": ")
                        } else {
                           pre <- ""
                        }
                        if ( !is.null(allNam[["group"]]) && !identical (allNam[["group"]], "wholeGroup") ) {
                           post <- zeile[allNam[["group"]]]
                           post <- paste(names(post), post, sep= "=", collapse=", ")
                        } else {
                           post <- ""
                        }
                        comps <- paste0(pre, post)
                     } else {
                        if(length(grep("groupdiff", zeile[["comparison"]], ignore.case=TRUE))==0){
                           IN <- ""
                        } else {                                                
                           if(length(grep("crossdiff", zeile[["comparison"]], ignore.case=TRUE))>0){
                              IN <- paste(gv, strsplit(zeile[[gv]], " - ")[[1]], collapse = " - ", sep="=")
                           } else {
                              IN <- eatTools::halveString(zeile[allNam[["group"]]], " - ")
                              for ( i in 1:nrow(IN)) {if(is.na(IN[i,2])) {IN[i,2]  <- IN[i,1]}}
                              IN <- t(IN)
                              IN <- paste(paste0("(", paste(colnames(IN), IN[1,], collapse=", ", sep="="), ")"), paste0("(", paste(colnames(IN), IN[2,], collapse=", ", sep="="), ")"), sep=" - ")
                           }
                        }
                        if(length(grep("crossdiff", zeile[["comparison"]], ignore.case=TRUE))>0){
                           cv    <- zeile[allNam[["group"]]]
                           weg   <- setdiff(grep(" - ", cv), grep("total", cv))
                           if ( length(weg)>0) {                                
                               cv    <- names(cv[-weg])                         
                           }  else  {
                               cv    <- names(cv)
                           }
                           cd    <- zeile[cv]
                           comps <- eatTools::makeDataFrame(eatTools::halveString(cd, " - "), verbose=FALSE)
                           for ( i in 1:nrow(comps)) {
                               if(is.na(comps[i,2])) {comps[i,2] <- comps[i,1]}
                               for ( j in 1:ncol(comps)) {comps[i,j] <- paste(rownames(comps)[i], comps[i,j], sep="=")} }
                           comps <- lapply(comps, FUN = function (co) {paste(co, collapse = ", ")})
                           if (IN != "") {
                               comps <- lapply(comps, FUN = function (co) {paste0("(", co, ": ",IN, ")")})
                           } else {
                               comps <- paste0("(", comps, ")")
                           }
                        } else {
                           comps <- IN
                        }
                        if(length(grep("trend", zeile[["comparison"]], ignore.case=TRUE))>0){
                           if( IN == "" && length(grep("crossdiff", zeile[["comparison"]], ignore.case=TRUE)) == 0) {
                               kl1 <- ""; kl2 <- ""; kl3 <- ""
                           } else {
                               kl1 <- "["; kl2 <- "]"; kl3 <- ": "
                           }
                           comps <- paste(paste0(paste0(kl1, paste(tv, strsplit(zeile[[tv]], " - ")[[1]], sep="="), kl3), comps, kl2), collapse=" - ")
                           if(zeile[["comparison"]] == "trend") {
                              if(length(allNam[["group"]])>0 && !identical (allNam[["group"]], "wholeGroup") ) {
                                 comps <- paste0(comps, ": ", paste(names(zeile[allNam[["group"]]]), zeile[allNam[["group"]]], sep="=", collapse = ", "))
                              } else {
                                 comps <- paste0(comps, ": total")
                              }
                           }
                        } else {
                           comps <- paste(comps, collapse=" - ")
                           if(!is.null(tv)) {
                              comps <- paste0("[", tv, "=",zeile[[tv]], ": ",comps, "]")
                           }
                        }
                     }
                     return(comps)})
           return(label2)}
           
report <- function ( repFunOut, trendDiffs = deprecated(), add=list(), exclude = c("NcasesValid", "var"), printGlm = FALSE,
                     round = TRUE, digits = 3, printDeviance = FALSE, printSE_correction = FALSE) {
     lifecycle::deprecate_warn("0.15.0", "report()", details = c(i = "For the original behavior of report() please use eatRep version 0.14.7: 'https://cran.r-project.org/src/contrib/Archive/eatRep/'"))
     if(!missing(trendDiffs)) {
         lifecycle::deprecate_warn("0.15.0", "report(trendDiffs)", details = c(i = "As differences of trends are equivalent to trends of differences, please look for 'trend_groupDiff' or 'trend_crossDiff' in the 'comparison' column."))
     }
     allN <- repFunOut[["allNam"]]
     out  <- eval(parse(text=paste0("report2(",paste(formalArgs(report2), formalArgs(report2), collapse=", ", sep="="), ")")))[["plain"]]
     gvars<- intersect(c(allN[["group"]], allN[["trend"]]), colnames(out))
     if ( length(gvars)>0) { for ( gv in gvars) {out[,gv] <- car::recode(gsub(" - ", ".vs.", gsub(" - total$", ".vs.wholeGroup",out[,gv])), "'total'=NA")}}
     avars<- intersect(c("parameter", "modus", "depVar", "comparison", "label1", names(add), allN[["group"]], allN[["trend"]], "es", "est", "p", "se"), colnames(out))
     if ("chiSquareTest" %in% out[,"parameter"]) {
         cat("Chi sqare test results cannot be transferred to old report() structure and will be ignored. Please use report2() instead.")
         out <- out[which(out[,"parameter"] != "chiSquareTest"),]
     }
     if(is.null(allN[["trend"]])) {
         outW <- out[,avars]
     }  else  {
         weg  <- sort(unique(strsplit(paste(gsub(".vs.", " ", out[,allN[["trend"]]]), collapse= " "), " ")[[1]]), decreasing=TRUE)
         weg1 <- paste("for year", weg)
         weg2 <- paste("(", lapply(combinat::combn(weg, 2, simplify=FALSE),FUN = paste, collapse=" - "), ")", sep="")
         for ( w in c(weg1, weg2)) {out[,"label1"] <- eatTools::removePattern(out[,"label1"], w)}
         mvars<- intersect(c("es", "est", "p", "se"), colnames(out))
         outW <- eatTools::makeDataFrame(tidyr::pivot_wider(out[,avars], names_from = allN[["trend"]], values_from = mvars), verbose=FALSE)
     }
     colnames(outW) <- car::recode(colnames(outW), "'label1'='group'")
     return(outW)}

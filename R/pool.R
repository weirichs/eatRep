pool.means <- function (m, se, na.rm = FALSE) {
     if(is.list(m))   {if(length(m) == 1)  { stopifnot(length(se)==1); m <- unlist(m); se <- unlist(se)}}
     if(!is.list(m))  {                                                         
        pooled <- mice::pool.scalar(Q=m, U=se^2)
        pooled <- data.frame ( m.pooled = pooled$qbar, se.pooled = sqrt(pooled$t), df = pooled$df, stringsAsFactors = FALSE)
     }  else  {                                                                 
        if(!all(unlist(lapply(m, length)) == unlist(lapply(se, length)) ) ) {   
            if(!all ( unlist(lapply(se,length)) == 0)) { 
               cat(paste("Some coefficients without standard error estimates. Standard errors won't be pooled.\n"))
            }
            se <- lapply(m, FUN = function ( x ) { rep(NA,length(x)) })
        }    
        M      <- length(m)
        N      <- length(m[[1]])
        Q.m    <- lapply(m,mean)                                                
        Qbar   <- mean(unlist(Q.m))                                             
        U      <- mean(unlist(lapply(se, FUN = function ( x ) { mean(x^2) } ))) 
        MS.b   <- N/(M-1) * sum((unlist(Q.m) - Qbar)^2)                         
        MS.omeg<- 1/(M*(N-1)) *  sum(unlist(lapply(m, FUN = function ( x ) { sum((x-mean(x))^2) })))
        varT   <- U+1/N*(1+1/M)*MS.b + (1-1/N)*MS.omeg                          
        dfN    <- 1 / ( (1/N*(1+1/M)*MS.b / varT)^2 * 1/(M-1) + (  ((1-1/N)*MS.omeg )/varT)^2 * 1/ (M * (N-1)) )
        pooled <- data.frame ( m.pooled = Qbar, se.pooled = sqrt(varT), df = dfN, stringsAsFactors = FALSE)
     }
     return(list(summary=pooled))}


pool.R2 <- function ( r2, N, quiet = FALSE ) {
           if(!is.list(r2)) {r2 <- list(r2)}
           if (missing(N))  {
               if(quiet == FALSE ) {cat("No sample size given. Will not compute standard error of pooled R squared.\n")}
               N <- lapply(r2, FUN = function (x) { rep ( 1000, length( x ) ) } )
               mis.N <- TRUE
           }
           if(!is.list(N))  {N  <- list(N)}
           if (!missing(N)) {
               stopifnot(length(N) == length(r2) )
               mis.N <- FALSE
               stopifnot( all ( sapply(N, length) == sapply(r2, length) ) )
           }
           Q.i     <- lapply(r2, FUN = function (x) {0.5*log( (1 + sqrt(x)) / (1-sqrt(x))  )})
           Q.i.err <- lapply(N,  FUN = function (n) {1 / (n-3)})
           untransformed <- pool.means(m = Q.i, se = Q.i.err)$summary[c("m.pooled","se.pooled")]
           transformed   <- as.data.frame ( ((exp(2*untransformed)-1) / (exp(2*untransformed)+1) )^2)
           if(mis.N) {return(transformed[1])} else {return(transformed)} }


jk2.pool <- function ( datLong, allNam, forceSingularityTreatment, modus ) {    
            retList <- do.call("rbind", by(data = datLong, INDICES = datLong[, c("group","parameter")], FUN = function ( u ) {
               comp <- table(u[,"comparison"], useNA="ifany")
               stopifnot(length(comp)==1)
               if(u[1,"parameter"] == "chiSquareTest") {                        
                  chi  <- by(u, INDICES = u[,c(allNam[["nest"]] )], FUN = function ( uN ) { uN[which(uN[,"coefficient"] == "chi2"),"value"]})
                  degFr<- by(u, INDICES = u[,c(allNam[["nest"]] )], FUN = function ( uN ) { uN[which(uN[,"coefficient"] == "df"),"value"]})
                  if ( length(table(degFr)) != 1 ) {
                       cat(paste0("Warning for '",u[1,"group"],"': degrees of freedom vary between imputations. Min: ",min(unlist(degFr)),"; Max: ", max(unlist(degFr)),". Chi-square test will be skipped.\n"))
                  }                                                             
                  degFr<- unique(unlist(degFr))[1]                              
                  pool <- miceadds::micombine.chisquare ( dk = unlist(chi), df=degFr, display = FALSE)
                  ret  <- data.frame ( group = names(table(u[,"group"])), depVar = allNam[["dependent"]], modus=modus, comparison = names(comp), parameter = names(table(u[,"parameter"])), coefficient = c("D2statistic","chi2Approx", "df1", "df2", "p", "pApprox"), value = pool[c("D", "chisq.approx", "df", "df2", "p", "p.approx")], u[1,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE, row.names=NULL)
               }  else  {
                  uM   <- by(u, INDICES = u[,c(allNam[["nest"]] )], FUN = function ( uN ) { uN[which(uN[,"coefficient"] == "est"),"value"]})
                  uSE  <- by(u, INDICES = u[,c(allNam[["nest"]] )], FUN = function ( uN ) { uN[which(uN[,"coefficient"] == "se"),"value"]})
                  if ( "es" %in% u[,"coefficient"] ) {
                        uES  <- by(u, INDICES = u[,c(allNam[["nest"]] )], FUN = function ( uN ) { uN[which(uN[,"coefficient"] == "es"),"value"]})
                        esP  <- mean(unlist(lapply(uES, mean)))
                  }
                  if(u[1,"parameter"] %in% c("R2", "R2nagel")) {                
                     getNvalid <- datLong[ intersect( intersect(  which(datLong[,"group"] == u[1,"group"]), which( datLong[,"parameter"] == "Nvalid")), which( datLong[,"coefficient"] == "est") ) ,]
                     if(nrow(getNvalid)==0) {
                        pooled    <- t(pool.R2(r2 = u[,"value"]))
                     }  else  {
                        if (forceSingularityTreatment == TRUE) {
                            if(nrow(getNvalid) != nrow(u) )  {
                               paste("Inconsistent number of sample size replications and/or R^2 estimates. Try workaround.\n")
                               u <- u[which(u[,"coefficient"] == "est"),]
                               stopifnot(nrow(getNvalid) == nrow(u))
                               pooled    <- t(pool.R2(r2 = u[,"value"], N = getNvalid[,"value"]))
                            }   else  {
                               pooled    <- t(pool.R2(r2 = u[,"value"], N = getNvalid[,"value"]))
                            }
                        }  else  {
                           pooled    <- t(pool.R2(r2 = u[,"value"], N = getNvalid[,"value"]))
                        }
                     }
                  } else {
                     pooled <- pool.means(m = uM, se = uSE)$summary[c("m.pooled","se.pooled")]
                  }
                  ret    <- data.frame ( group = names(table(u[,"group"])), depVar = allNam[["dependent"]], modus=modus, comparison = names(comp), parameter = names(table(u[,"parameter"])), coefficient = c("est","se"), value = unlist(pooled), u[1,allNam[["group"]],drop=FALSE], stringsAsFactors = FALSE, row.names=NULL)
                  if ( "es" %in% u[,"coefficient"] ) {
                        retI <- ret[1,]
                        retI[,"coefficient"] <- "es"
                        retI[,"value"]       <- esP
                        ret  <- rbind(ret, retI)
                  }
               }
               return(ret)}))
           return(retList)}

pool.corr <- function( corrs , N , conf.level = .05){
        fisher.corrs <- 1/2*log( ( 1 + corrs) / ( 1 - corrs ) )                 
        var.fisher <- rep( 1/(N-3) , length(corrs) )                            
        fisher.cor.combine <- mice::pool.scalar( fisher.corrs , var.fisher)
        zr <- fisher.cor.combine$qbar
        zr.se <- sqrt( fisher.cor.combine$t )
        t.zr <- zr / zr.se
        fisher2cor <- function(z){ ( exp(2*z) - 1 )/ ( exp(2*z) + 1 ) }
        res <- c( "r" = fisher2cor(zr)  ,
            "fisher_r" = zr ,
            "fisher_rse" = zr.se ,
            "t" = t.zr  ,
            "p" = 2 * pnorm( abs(t.zr) , lower.tail = FALSE ) ,
             fisher2cor( zr + qnorm( ( 1 - conf.level ) / 2 ) * zr.se ) ,
             fisher2cor( zr - qnorm( ( 1 - conf.level ) / 2 ) * zr.se ) )
            names(res)[6] <- paste( "lower" , round(100*conf.level,2),sep="")
            names(res)[7] <- paste( "upper" , round(100*conf.level,2),sep="")
        res <- c( res , ( res[6] - res[7] ) / ( 2* qnorm( ( 1 - conf.level )/2 ) ) )
        names(res)[8] <- "rse"
        res <- res[ c(1,8,2:7) ]
        res <- round(res, 6)
        return(res) }           

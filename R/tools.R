make.indikator <- function(variable, name.var = "ind", force.indicators = NULL, separate.missing.indikator = c("no","ifany", "always"), sep = "_" )  {
                  separate.missing.indikator <- match.arg(separate.missing.indikator)
                  t.var <- table(variable, useNA = separate.missing.indikator )
                  if(!is.null(force.indicators))  {
                     additional <- setdiff(as.character(force.indicators), names(t.var))
                     if(length(additional) > 0 )  {
                        add <- rep(0,length(additional))
                        names(add) <- additional
                        t.var      <- c(t.var, add)
                     }
                  }
                  ind.i <- data.frame( variable, lapply(names(t.var), FUN = function(iii) {
                           if(!is.na(iii)) {ind.iii  <- which(variable == iii)}
                           if(is.na(iii))  {ind.iii  <- which(is.na(variable))}
                           ret      <- rep(0, length(variable) )
                           if(length(ind.iii)>0)  {ret[ind.iii] <- 1}
                           if(separate.missing.indikator == "no" )  {
                              if(length(which(is.na(variable))) > 0 ) {
                                 ret[which(is.na(variable))] <- NA
                              }
                           }
                           return(ret)}), stringsAsFactors = FALSE )
                  colnames(ind.i)[-1] <- paste(name.var, names(t.var), sep=sep)
                  return(ind.i)}

identifyMode <- function ( name, type) {
            res <- paste0(car::recode(type, "'NONE'='CONV'"),".", name)
            return(res)}

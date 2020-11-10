jk2.mean <- function ( ... ) {
     message("Function 'jk2.mean' was deprecated. Please use 'comp.stats' instead.")
     ret <- comp.stats ( ... )
     return(ret)
}
jk2.table <- function ( ... ) {
     message("Function 'jk2.table' was deprecated. Please use 'comp.table' instead.")
     ret <- comp.table ( ... )
     return(ret)
}
jk2.quantile <- function ( ... ) {
     message("Function 'jk2.quantile' was deprecated. Please use 'comp.quantile' instead.")
     ret <- comp.quantile ( ... )
     return(ret)
}
jk2.glm <- function ( ... ) {
     message("Function 'jk2.glm' was deprecated. Please use 'comp.glm' instead.")
     ret <- comp.glm ( ... )
     return(ret)
}


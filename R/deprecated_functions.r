jk2.mean <- function ( ... ) {
     message("Function 'jk2.mean' was deprecated. Please use 'repMean' instead.")
     ret <- repMean ( ... )
     return(ret)
}
jk2.table <- function ( ... ) {
     message("Function 'jk2.table' was deprecated. Please use 'repTable' instead.")
     ret <- repTable ( ... )
     return(ret)
}
jk2.quantile <- function ( ... ) {
     message("Function 'jk2.quantile' was deprecated. Please use 'repQuantile' instead.")
     ret <- repQuantile ( ... )
     return(ret)
}
jk2.glm <- function ( ... ) {
     message("Function 'jk2.glm' was deprecated. Please use 'repGlm' instead.")
     ret <- repGlm ( ... )
     return(ret)
}


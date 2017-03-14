ranef <- function(x, ...) UseMethod("ranef")

ranef.mexhaz <- function(x, ...){

    if (is.na(x$random)){
        stop("The ranef.mexhaz() function can only be used if a random effect was included in the model...")
    }
    se.mu.hat <- sqrt(diag(x$var.mu.hat))

    TAB <- cbind(x$mu.hat,se.mu.hat)

    return(TAB)
}

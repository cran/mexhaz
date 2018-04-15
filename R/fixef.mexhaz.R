fixef <- function(x, ...) UseMethod("fixef")

fixef.mexhaz <- function(x, ...){

    if (is.na(x$random)){
        res <- x$coefficients[1:x$n.par]
    }
    else {
        res <- x$coefficients[1:(x$n.par-1)]
    }

    return(res)
}

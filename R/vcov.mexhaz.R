vcov <- function(object, ...) UseMethod("vcov")

vcov.mexhaz <- function(object, ...){

    if (is.na(object$random)){
        res <- object$vcov
    }
    else {
        res <- object$vcov[1:(object$n.par-1),1:(object$n.par-1)]
    }

    return(res)
}

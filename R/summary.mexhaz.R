summary.mexhaz <- function(object, ...){
    coef <- coef(object)
    se <- sqrt(diag(object$vcov))
    tval <- coef(object)/se
    df <- object$n.obs-object$n.par
    n.miss <- object$n.obs.tot-object$n.obs
    IdxPH <- which(names(coef)%in%object$names.ph)

    TAB <- cbind(Estimate = coef,
                 StdErr = se,
                 t.value = tval,
                 p.value = 2*pt(-abs(tval),df=df))

    x <- list(call=object$call,
              coefficients=TAB,
              n.obs=object$n.obs,
              n.events=object$n.events,
              n.miss=(object$n.obs.tot-object$n.obs),
              n.time.0=object$n.time.0,
              n.par=object$n.par,
              loglik=object$loglik,
              df=object$n.obs-object$n.par,
              idx.ph=IdxPH)
    class(x) <- "summary.mexhaz"
    x
}

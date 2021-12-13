summary.mexhaz <- function(object, ...){
    coef <- coef(object)
    se <- sqrt(diag(object$vcov))
    zval <- coef(object)/se
    n.miss <- object$n.obs.tot-object$n.obs
    IdxPH <- which(names(coef)%in%object$names.ph)

    TAB <- cbind(Estimate = coef,
                 StdErr = se,
                 z.value = zval,
                 p.value = 2*(1-pnorm(abs(zval))))

    x <- list(call=object$call,
              coefficients=TAB,
              n.obs=object$n.obs,
              n.events=object$n.events,
              n.miss=(object$n.obs.tot-object$n.obs),
              n.time.0=object$n.time.0,
              n.par=object$n.par,
              loglik=object$loglik,
              idx.ph=IdxPH)
    class(x) <- "summary.mexhaz"
    x
}

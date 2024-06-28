predict.mexhaz <- function (object, time.pts, data.val = data.frame(.NotUsed = NA),
    marginal = FALSE, quant.rdm = 0.5, cluster = NULL, conf.int = c("delta", "simul", "none"), level = 0.95,
    delta.type.h = c("log", "plain"), delta.type.s = c("log-log", "log", "plain"), nb.sim = 10000,
    keep.sim = FALSE, include.gradient = FALSE, dataset = NULL, ...)
{
    if (object$recurrent==FALSE){
        res <- predictStd(object=object, time.pts=time.pts, data.val=data.val,
                          marginal=marginal, quant.rdm=quant.rdm, cluster=cluster, conf.int=conf.int, level=level,
                          delta.type.h=delta.type.h, delta.type.s=delta.type.s, nb.sim=nb.sim, keep.sim=keep.sim,
                          include.gradient=include.gradient, dataset=dataset, ...)
    }
    else {
        res <- predictRec(object=object, time.pts=time.pts, data.val=data.val,
                          conf.int=conf.int, level=level, nb.sim=nb.sim, keep.sim=keep.sim, include.gradient=include.gradient, ...)
    }
    return(res)
}

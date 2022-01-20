points.predMexhaz <- function(x,which=c("surv","hazard"),conf.int=TRUE,lty.pe="solid",lty.ci="dashed",...){
    which <- match.arg(which)
    if (x$type=="multiobs"){
        stop("The 'points.predMexhaz' function applies only to predictions realised on a single vector of covariates.")
    }
    time.pts <- x$results$time.pts
    if (which=="hazard"){
        points(time.pts,x$results$hazard,lty=lty.pe,...)
    }
    else {
        points(c(0,time.pts),c(1,x$results$surv),lty=lty.pe,...)
    }
    if (conf.int==TRUE & x$ci.method!="none"){
        if (which=="hazard"){
            points(time.pts,x$results$hazard.inf,lty=lty.ci,...)
            points(time.pts,x$results$hazard.sup,lty=lty.ci,...)
        }
        else {
            points(c(0,time.pts),c(1,x$results$surv.inf),lty=lty.ci,...)
            points(c(0,time.pts),c(1,x$results$surv.sup),lty=lty.ci,...)
        }
    }
}

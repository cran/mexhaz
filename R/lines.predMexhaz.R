lines.predMexhaz <- function(x,which=c("surv","hazard"),conf.int=TRUE,lty.pe="solid",lty.ci="dashed",...){
    which <- match.arg(which)
    if (x$type=="multiobs"){
        stop("The 'lines.predMexhaz' function applies only to predictions realised on a single vector of covariates.")
    }
    time.pts <- x$results$time.pts
    if (which=="hazard"){
        lines(time.pts,x$results$hazard,lty=lty.pe,...)
    }
    else {
        lines(c(0,time.pts),c(1,x$results$surv),lty=lty.pe,...)
    }
    if (conf.int==TRUE & x$ci.method!="none"){
        if (which=="hazard"){
            lines(time.pts,x$results$hazard.inf,lty=lty.ci,...)
            lines(time.pts,x$results$hazard.sup,lty=lty.ci,...)
        }
        else {
            lines(c(0,time.pts),c(1,x$results$surv.inf),lty=lty.ci,...)
            lines(c(0,time.pts),c(1,x$results$surv.sup),lty=lty.ci,...)
        }
    }
}

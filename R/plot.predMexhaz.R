plot.predMexhaz <- function(x,which=c("surv","hazard"),conf.int=TRUE,lty.pe="solid",lty.ci="dashed",...){
    which <- match.arg(which)
    if (x$type=="multiobs"){
        stop("The 'plot.predMexhaz' function applies only to predictions realised on a single vector of covariates.")
    }
    time.pts <- x$results$time.pts
    if (which=="hazard"){
        plot(time.pts,x$results$hazard,type="l",xaxs="i",xlab="Time",ylab="Hazard",lty=lty.pe,...)
    }
    else {
        plot(c(0,time.pts),c(1,x$results$surv),type="l",xaxs="i",xlab="Time",ylab="Survival",lty=lty.pe,...)
    }
    if (conf.int==TRUE & x$ci.method!="none"){
        if (which=="hazard"){
            points(time.pts,x$results$hazard.inf,type="l",lty=lty.ci,...)
            points(time.pts,x$results$hazard.sup,type="l",lty=lty.ci,...)
        }
        else {
            points(c(0,time.pts),c(1,x$results$surv.inf),type="l",lty=lty.ci,...)
            points(c(0,time.pts),c(1,x$results$surv.sup),type="l",lty=lty.ci,...)
        }
    }
}

lines.resMexhaz <- function(x,conf.int=TRUE,lty.pe="solid",lty.ci="blank",col.ci="blue",alpha.col.ci=0.25,...){
    lab <- c("Hazard Ratio","Risk Ratio","Survival Difference","Conditional Survival","Direct Adjusted Survival")[which(c("hr","rr","sd","cs","as")%in%x$type)]
    if (x$multiobs){
        stop("The 'lines.resMexhaz' function applies only to predictions realised on a single vector of covariates.")
    }
    time.pts <- x$results$time.pts
    lines(time.pts,x$results[,2],lty=lty.pe,...)
    if (conf.int==TRUE){
        polygon(c(time.pts,rev(time.pts)),c(x$results[,3],rev(x$results[,4])),col=rgb(t(col2rgb(col.ci))/255,alpha=alpha.col.ci),border=NA)
        lines(time.pts,x$results[,2],lty=lty.pe,...)
        lines(time.pts,x$results[,3],lty=lty.ci,...)
        lines(time.pts,x$results[,4],lty=lty.ci,...)
    }
}

print.summary.mexhaz <- function(x, ...){
    HRTab <- NULL
    ishr <- 0
    if (length(x$idx.ph)>0){
        PH <- x$coefficients[x$idx.ph,1]
        sePH <- x$coefficients[x$idx.ph,2]
        HR <- exp(PH)
        ci.lower <- exp(PH+qt(0.025,df=x$df)*sePH)
        ci.upper <- exp(PH+qt(0.975,df=x$df)*sePH)
        HRTab <- cbind(Coef=round(PH,4),HR=round(HR,4),CI.lower=round(ci.lower,4),CI.upper=round(ci.upper,4))
        ishr <- 1
    }

    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
    if (ishr){
        cat("\nHazard ratios (for proportional effect variables):\n")
        print(HRTab)
    }
    if(x$n.miss>0){
        Miss <- paste(x$n.miss," observations were deleted due to missingness\n",sep="")
    }
    else {
        Miss <- NULL
    }
    if (x$n.time.0>0){
        Time0 <- paste(x$n.time.0," observations had a follow-up time equal to 0 (replaced by 1/730.5)\n",sep="")
    }
    else {
        Time0 <- NULL
    }

    cat(paste("\nlog-likelihood: ",round(x$loglik,4)," (for ",x$n.par," degree(s) of freedom)\n",sep=""))
    cat(paste("\nnumber of observations: ",x$n.obs,", number of events: ",x$n.events,"\n",sep=""))
    cat(Miss)
    cat(Time0)
}

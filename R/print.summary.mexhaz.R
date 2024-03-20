print.summary.mexhaz <- function(x, ...){
    HRTab <- NULL
    ishr <- 0
    if (length(x$idx.ph)>0){
        PH <- x$coefficients[x$idx.ph,1]
        sePH <- x$coefficients[x$idx.ph,2]
        HR <- exp(PH)
        ci.lower <- exp(PH+qnorm(0.025)*sePH)
        ci.upper <- exp(PH+qnorm(0.975)*sePH)
        HRTab <- cbind(Coef=round(PH,4),HR=round(HR,4),CI.lower=round(ci.lower,4),CI.upper=round(ci.upper,4))
        row.names(HRTab) <- x$names.ph
        ishr <- 1
    }

    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
    if (x$baseline=="pw.cst"){
        cat("\nNote: For a piecewise constant baseline hazard, the coefficients\n correspond to the logarithm of the hazard on each time interval.\n")
    }
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

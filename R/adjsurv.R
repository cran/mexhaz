adjsurv <- function(object,time.pts,data,data.0=NULL,weights=NULL,marginal=TRUE,quant.rdm=0.5,cluster=NULL,quant.rdm.0=0.5,cluster.0=NULL,level=0.95,dataset=NULL){

    if (is.na(object$random)){
        marginal <- FALSE
    }
    Idx.T.NA <- which(is.na(time.pts) | time.pts<0)
    if (length(Idx.T.NA)>0){
        stop("The 'time.pts' argument contains NA or negative values...")
    }
    Idx.T.Max <- which(time.pts>object$max.time)
    if (length(Idx.T.Max)>0){
        warning(paste("The model cannot be used to predict survival for times greater than ",object$max.time," (maximum follow-up time on which the object estimation was based). Consequently, these time values were removed from the 'time.pts' vector.",sep=""))
        time.pts <- time.pts[-Idx.T.Max]
    }
    which.zero <- which(time.pts==0)
    time.pts <- time.pts[time.pts!=0]
    nb.time.pts <- length(time.pts)
    if (nb.time.pts==0){
        stop("The 'time.pts' argument contains no values for which predictions can be made...")
    }
    if (!is.numeric(level) | (level>1 | level<0)){
        level <- 0.95
        warning("The 'level' argument must be a numerical value in (0,1)...")
    }
    alpha <- (1-level)/2

    NbObs <- dim(data)[1]
    if (is.null(weights)){
        weights <- rep(1/NbObs,NbObs)
    }
    if (length(weights)!=NbObs){
        stop("The argument 'weights' must be NULL or a vector of length the number of rows in 'data'...")
    }
    if (!is.null(data.0)){
        if (dim(data.0)[1]!=NbObs){
            stop("'data.0' must contain the same number of observations as 'data'...")
        }
        if (is.null(cluster.0)){
            cluster.0 <- cluster
        }
    }

    SurvPop <- VarSurv <- rep(0,nb.time.pts)
    if (!is.null(data.0)){
        DiffSurv <- SurvPop0 <- VarDiff <- VarSurv0 <- rep(0,nb.time.pts)
    }

    for (i in 1:nb.time.pts){
        PredI <- predict.mexhaz(object,time.pts[i],data.val=data,cluster=cluster,marginal=marginal,quant.rdm=quant.rdm,include.gradient=TRUE,dataset=dataset)
        surv <- PredI$results$surv
        SurvPop[i] <- weights%*%surv
        wgt.grad.surv <- -t(PredI$grad.logcum*surv*log(surv))%*%weights
        VarSurv[i] <- t(wgt.grad.surv)%*%PredI$vcov%*%wgt.grad.surv
        if (!is.null(data.0)){
            PredI0 <- predict.mexhaz(object,time.pts[i],data.val=data.0,cluster=cluster.0,marginal=marginal,quant.rdm=quant.rdm.0,include.gradient=TRUE,dataset=dataset)
            surv0 <- PredI0$results$surv
            SurvPop0[i] <- weights%*%surv0
            wgt.grad.surv0 <- -t(PredI0$grad.logcum*surv0*log(surv0))%*%weights
            VarSurv0[i] <- t(wgt.grad.surv0)%*%PredI0$vcov%*%wgt.grad.surv0
            DiffSurv[i] <- weights%*%(surv-surv0)
            wgt.grad.diffsurv <- -t(PredI$grad.logcum*surv*log(surv)-PredI0$grad.logcum*surv0*log(surv0))%*%weights
            VarDiff[i] <- t(wgt.grad.diffsurv)%*%PredI$vcov%*%wgt.grad.diffsurv
        }
    }
    LogCumPop <- log(-log(SurvPop))
    VarLCP <- VarSurv/(SurvPop*log(SurvPop))^2
    SP <- exp(-exp(LogCumPop+sqrt(VarLCP)%*%t(c(0,qnorm(alpha),qnorm(1-alpha)))))
    if (length(which.zero)>0){
        time.pts <- c(0,time.pts)
        SP <- rbind(1,SP)
    }
    if (!is.null(data.0)){
        LogCumPop0 <- log(-log(SurvPop0))
        VarLCP0 <- VarSurv0/(SurvPop0*log(SurvPop0))^2
        SP0 <- exp(-exp(LogCumPop0+sqrt(VarLCP0)%*%t(c(0,qnorm(alpha),qnorm(1-alpha)))))
        DSP <- DiffSurv+sqrt(VarDiff)%*%t(c(0,qnorm(alpha),qnorm(1-alpha)))
        if (length(which.zero)>0){
            SP0 <- rbind(1,SP0)
            DSP <- rbind(0,DSP)
        }
        results <- as.data.frame(cbind(time.pts,SP,SP0,DSP))
        names(results) <- c("time.pts","adj.surv","adj.surv.inf","adj.surv.sup","adj.surv.0","adj.surv.0.inf","adj.surv.0.sup","diff.adj.surv","diff.adj.surv.inf","diff.adj.surv.sup")
    }
    else {
        results <- as.data.frame(cbind(time.pts,SP))
        names(results) <- c("time.pts","adj.surv","adj.surv.inf","adj.surv.sup")
    }
    res.AS <- list(results=results,type="as",multiobs=FALSE,ci.method="delta",level=level)
    class(res.AS) <- "resMexhaz"
    return(res.AS)

}

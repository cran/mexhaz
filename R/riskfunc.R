riskfunc <- function(object,time.pts,data,data.0,clust.name=NULL,marginal=TRUE,type=c("hr","rr"),conf.int=c("delta","simul"),level=0.95,nb.sim=10000,seed=NULL){
    if (is.null(seed)){
        seed <- .Random.seed[1]
    }
    if (is.null(time.pts)){
        time.pts <- object$max.time/2
    }

    dim.data <- dim(data)[1]
    if (dim(data.0)[1]!=dim.data){
        stop("Arguments 'data' and 'data.0' must have the same number of rows...")
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
    if (length(time.pts)==0){
        stop("The 'time.pts' argument contains no values for which predictions can be made...")
    }
    nb.time.pts <- length(time.pts)
    if (nb.time.pts>1 & dim.data>1) {
        stop("Predictions can be made for n individuals at 1 time point or for 1 individual at m time points but not for several individuals at several time points...")
    }
    multiobs <- (dim.data>1)
    cluster.0 <- cluster.1 <- NULL
    if (!is.null(clust.name)){
        if (!is.character(clust.name[1]) & length(clust.name)!=1){
            stop("The 'clust.name' argument must be a character string giving the name of the clustering variable in the dataset 'data'/'data.0'...")
        }
        if (!clust.name%in%names(data)){
            stop("The 'clust.name' argument must be a character string giving the name of the clustering variable in the dataset 'data'/'data.0'...")
        }
        cluster.1 <- data[,clust.name]
        cluster.0 <- data.0[,clust.name]
    }
    if (!is.numeric(level) | (level>1 | level<0)){
        level <- 0.95
        warning("The 'level' argument must be a numerical value in (0,1)...")
    }
    alpha <- (1-level)/2

    keep.sim <- (conf.int=="simul")

    set.seed(seed)
    p0 <- predict.mexhaz(object,time.pts=time.pts,data.val=data.0,cluster=cluster.0,marginal=marginal,conf.int=conf.int,include.gradient=TRUE,nb.sim=nb.sim,keep.sim=keep.sim)
    set.seed(seed)
    p1 <- predict.mexhaz(object,time.pts=time.pts,data.val=data,cluster=cluster.1,marginal=marginal,conf.int=conf.int,include.gradient=TRUE,nb.sim=nb.sim,keep.sim=keep.sim)
    time.pts <- p0$results$time.pts
    quant <- t(c(0,qnorm(alpha),qnorm(1-alpha)))

    if (type=="hr"){
        LogHR <- log(p1$results$hazard) - log(p0$results$hazard)
        if (conf.int=="delta"){
            Grad.LogHR <- p1$grad.loghaz - p0$grad.loghaz
            Var.LogHR <- diag(Grad.LogHR%*%p1$vcov%*%t(Grad.LogHR))
            HR <- exp(LogHR+sqrt(Var.LogHR)%*%quant)
        }
        else if (conf.int=="simul"){
            sim.HR <- exp(log(p1$sim.haz)-log(p0$sim.haz))
            BInf.HR <- apply(sim.HR,1,quantile,prob=alpha)
            BSup.HR <- apply(sim.HR,1,quantile,prob=1-alpha)
            HR <- cbind(exp(LogHR),BInf.HR,BSup.HR)
        }
        results <- as.data.frame(cbind(time.pts,HR))
        names(results) <- c("time.pts","hazard.ratio","hr.inf","hr.sup")
    }
    else if (type=="rr"){
        S1 <- p1$results$surv
        S0 <- p0$results$surv
        LogRR <- log(1-S1)-log(1-S0)
        if (conf.int=="delta"){
            Grad.LogRR <- (S0*log(S0)/(1-S0))*p0$grad.logcum - (S1*log(S1)/(1-S1))*p1$grad.logcum
            Var.LogRR <- diag(Grad.LogRR%*%p1$vcov%*%t(Grad.LogRR))
            RR <- exp(LogRR+sqrt(Var.LogRR)%*%quant)
        }
        else if (conf.int=="simul"){
            sim.S1 <- p1$sim.surv
            sim.S0 <- p0$sim.surv
            sim.RR <- exp(log(1-sim.S1)-log(1-sim.S0))
            BInf.RR <- apply(sim.RR,1,quantile,prob=alpha,na.rm=TRUE)
            BSup.RR <- apply(sim.RR,1,quantile,prob=1-alpha,na.rm=TRUE)
            RR <- cbind(exp(LogRR),BInf.RR,BSup.RR)
        }
        results <- as.data.frame(cbind(time.pts,RR))
        names(results) <- c("time.pts","risk.ratio","rr.inf","rr.sup")
    }
    res.RR <- list(results=results,type=type,multiobs=multiobs,ci.method=conf.int,level=level)
    class(res.RR) <- "resMexhaz"
    return(res.RR)
}

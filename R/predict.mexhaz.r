predict.mexhaz <- function(object,time.pts,data.val=data.frame(.NotUsed=NA),cluster=NULL,conf.int=c("delta","simul","none"),level=0.95,delta.type.h=c("log","plain"),delta.type.s=c("log-log","log","plain"),nb.sim=10000,include.gradient=FALSE,...){

    time0 <- as.numeric(proc.time()[3])
    conf.int <- match.arg(conf.int)
    delta.type.h <- match.arg(delta.type.h)
    delta.type.s <- match.arg(delta.type.s)
    call <- object$call
    formula <- object$formula
    Xlev <- object$xlevels
    base <- object$base
    degree <- object$degree
    knots <- object$knots
    n.gleg <- object$n.gleg
    BoI <- object$boundary.knots[1]
    BoS <- object$boundary.knots[2]
    max.time <- object$max.time
    coef <- object$coefficients
    vcov <- object$vcov
    n.par <- object$n.par
    if (!is.numeric(level) | (level>1 | level<0)){
        level <- 0.95
        warning("The 'level' argument must be a numerical value in (0,1)...")
    }
    if (!is.null(cluster)){
        if (is.na(object$random)){
            stop("The 'cluster' argument cannot be used with a fixed effects model...")
        }
        if (length(cluster)!=1){
            stop("The 'cluster' argument must be a single number or character string corresponding to the name of the cluster for which predictions are wanted...")
        }
        Idx.mh <- which(object$mu.hat$cluster==cluster)
        if (length(Idx.mh)!=1){
            stop(paste("There is no cluster named '",cluster,"'...",sep=""))
        }
        coef <- c(coef[-n.par],object$mu.hat[Idx.mh,]$mu.hat)
        vcov1 <- cbind(vcov[-n.par,-n.par],t(object$vcov.fix.mu.hat[Idx.mh,,drop=FALSE]))
        vcov <- rbind(vcov1,c(object$vcov.fix.mu.hat[Idx.mh,],object$var.mu.hat[Idx.mh,Idx.mh]))
    }
    Idx.T.NA <- which(is.na(time.pts) | time.pts<0)
    if (length(Idx.T.NA)>0){
        stop("The 'time.pts' argument contains NA or negative values...")
    }
    Idx.T.Max <- which(time.pts>max.time)
    if (length(Idx.T.Max)>0){
        warning(paste("The model cannot be used to predict survival for times greater than ",max.time," (maximum follow-up time on which the model estimation was based). Consequently, these time values were removed from the 'time.pts' vector.",sep=""))
        time.pts <- time.pts[-Idx.T.Max]
    }
    time.pts <- time.pts[time.pts!=0]
    if (length(time.pts)==0){
        stop("The 'time.pts' argument contains no values for which predictions can be made...")
    }
    time.pts <- time.pts[order(time.pts)]
    nb.time.pts <- length(time.pts)
    if (nb.time.pts>1 & dim(data.val)[1]>1) {
        stop("Predictions can be made for n individuals at 1 time point or for 1 individual at m time points but not for several individuals at several time points...")
    }

    if (nb.time.pts>1){
        typepred <- "multitime"
    }
    else {
        typepred <- "multiobs"
    }

    if (conf.int=="simul" & (nb.sim<=0 | (round(nb.sim,0)!=nb.sim))){
        warning("The 'nb.sim' argument must be a strictly positive integer. Consequently, the default value of 50000 simulations has been used...")
        nb.sim <- 50000
    }

    tot.formula <- terms(formula,data=data.val,specials="nph")
    indNph <- attr(tot.formula,"specials")$nph
    if (length(indNph)>0){
        nTerm <- NULL
        nTerm2 <- NULL
        for (i in 1:length(indNph)){
            nphterm <- attr(tot.formula,"variables")[[1+indNph[i]]]
            nTerm <- c(nTerm,deparse(nphterm[[2L]],width.cutoff=500L,backtick=TRUE))
            nTerm2 <- c(nTerm2,deparse(nphterm,width.cutoff=500L,backtick=TRUE))
        }
        nTerm <- paste(nTerm,collapse="+")
        nTerm2 <- paste(nTerm2,collapse="-")
        FormulaN <- as.formula(paste("~",nTerm))
    }
    else {
        nTerm2 <- 0
        FormulaN <- as.formula("~1")
    }
    data.new.temp <- cbind(time.pts, data.val, row.names = NULL)
    FormulaF <- update(tot.formula,paste("NULL~.",nTerm2,sep="-"))
    Test <- model.frame(FormulaF,data=data.new.temp,xlev=Xlev,na.action=NULL) # Will raise an error if one of the variables necessary to the model is missing from data.val
    if (dim(Test)[1]==0){
        stop("The formula is incompatible with the data provided. If bs() was used to model the effect of some of the variables included in the formula, make sure that the spline boundaries were specified...")
    }

    Idx.Dv.NA <- which(apply(Test,1,function(x){sum(is.na(x))})>0)
    if (length(Idx.Dv.NA)>0){
        warning("Some rows in the 'data.val' data.frame had NA values. Consequently, these rows were deleted.")
        Test <- Test[-Idx.Dv.NA,,drop=FALSE]
    }
    if (dim(Test)[1]==0){
        stop("The 'data.val' dataframe contains no valid row on which predictions can be based.")
    }
    data.new <- cbind(time.pts,Test,row.names=NULL)
    time.pts <- data.new[,1]
    nb.time.pts <- length(time.pts)
    FALCenv <- environment()

    if (base=="pw.cst"){
        HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
            .Call(C_HazardBs0R,x,nph,timecat,fixobs,param,paramf,as.double(matk))
        }
        DeltaFct <- function(x0,x,nph,timecat0,timecat,fixobs,paramt,deg,n,lw,matk,totk,intk,nsadj1,nsadj2,varcov,grad){
            .Call(C_DeltaBs0R,x,nph,timecat,fixobs,paramt,as.double(matk),as.double(varcov),grad)
        }
    }
    else if (base=="exp.bs" & degree==1){
        HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
            .Call(C_HazardBs1R,x,nph,timecat,fixobs,param,paramf,as.double(matk),as.double(totk))
        }
        DeltaFct <- function(x0,x,nph,timecat0,timecat,fixobs,paramt,deg,n,lw,matk,totk,intk,nsadj1,nsadj2,varcov,grad){
            .Call(C_DeltaBs1R,x,nph,timecat,fixobs,paramt,as.double(matk),as.double(totk),as.double(varcov),grad)
        }
    }
    else if (base=="exp.bs" & degree%in%c(2,3)){
        HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
            .Call(C_HazardBs23R,x,nph,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
        }
        DeltaFct <- function(x0,x,nph,timecat0,timecat,fixobs,paramt,deg,n,lw,matk,totk,intk,nsadj1,nsadj2,varcov,grad){
            .Call(C_DeltaBs23R,x,nph,timecat,fixobs,paramt,deg,n,lw,as.double(matk),as.double(totk),as.double(varcov),grad)
        }
    }
    else if (base=="exp.ns"){
        HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
            .Call(C_HazardNsR,x,nph,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2))
        }
        DeltaFct <- function(x0,x,nph,timecat0,timecat,fixobs,paramt,deg,n,lw,matk,totk,intk,nsadj1,nsadj2,varcov,grad){
            .Call(C_DeltaNsR,x,nph,timecat,fixobs,paramt,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2),as.double(varcov),grad)
        }
    }
    else if (base=="weibull"){
        HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
            log.p.LT.1 <- log(paramf[1]) + as.vector(fixobs%*%paramf[-1])
            log.p.LT.2 <- log(param[1]) + as.vector(nph%*%param[-1])
            l.lambda.beta <- log.p.LT.2 +x*(exp(log.p.LT.2)-1)+log.p.LT.1
            l.Lambda.beta <- x*exp(log.p.LT.2)+log.p.LT.1
            Result <- list(LogHaz=l.lambda.beta, LogCum=l.Lambda.beta)
            return(Result)
        }
        DeltaFct <- function(x0,x,nph,timecat0,timecat,fixobs,paramt,deg,n,lw,matk,totk,intk,nsadj1,nsadj2,varcov,grad){
            .Call(C_DeltaWeibR,x,t(nph),t(fixobs),paramt,varcov,grad)
        }
    }
    Surv <- NA
    lambda <- NA
    varcov <- NA
    BInf1 <- NA
    BInf2 <- NA
    BSup1 <- NA
    BSup2 <- NA
    Var.Log.Haz <- NA
    Var.Log.Cum <- NA
    Grad.LH <- NA
    Grad.LC <- NA


    # Functions for computing the part of the B-spline bases that depends only on the knots and can therefore be calculated only once at the beginning of the function (used in combination with the IntHazard function to estimate the hazard and the cumulative hazard)

    # For B-splines (also used for Natural Cubic Splines)
    TransfDeg <- function(vec.knots,deg){
        if (deg==1){
            Le <- length(vec.knots)
            Res <- vec.knots[2:Le]-vec.knots[1:(Le-1)]
        }
        else {
            Dim <- (length(vec.knots)-(2*deg-1))
            Res <- matrix(NA,2*(deg-1),Dim)
            if (deg==2){
                for (i in 1:Dim){
                    TempK <- vec.knots[i:(i+3)]
                    Res[1,i] <- 1/((TempK[4]-TempK[2])*(TempK[3]-TempK[2]))
                    Res[2,i] <- 1/((TempK[3]-TempK[1])*(TempK[3]-TempK[2]))
                }
            }
            else if (deg==3){
                for (i in 1:Dim){
                    TempK <- vec.knots[i:(i+5)]
                    Res[1,i] <- 1/((TempK[6]-TempK[3])*(TempK[5]-TempK[3])*(TempK[4]-TempK[3]))
                    Res[2,i] <- 1/((TempK[5]-TempK[2])*(TempK[4]-TempK[2])*(TempK[4]-TempK[3]))
                    Res[3,i] <- 1/((TempK[5]-TempK[2])*(TempK[5]-TempK[3])*(TempK[4]-TempK[3]))
                    Res[4,i] <- 1/((TempK[4]-TempK[1])*(TempK[4]-TempK[2])*(TempK[4]-TempK[3]))
                }
            }
        }
        return(Res)
    }

    # For Natural Cubic Splines

    if (base=="exp.bs"){
        NsAdjust <- function(vec.knots,BoI,BoS){
            NULL
        }
        dbase <- degree
    }
    else if (base=="exp.ns"){
        NsAdjust <- function(vec.knots,BoI,BoS){
            Dim <- (length(vec.knots)-5)
            Diag <- diag(Dim)
            QR <- qr(t(splineDesign(vec.knots,c(BoI,BoS),4,c(2,2))[,-1,drop=FALSE]))
            Res1 <- t(apply(Diag,1,function(x){qr.qty(QR,x)})[-c(1:2),,drop=FALSE])
            SpI <- c(BoI,splineDesign(vec.knots,rep(BoI,2L),4,c(0,1))[2,1:2])
            SpS <- c(BoS,splineDesign(vec.knots,rep(BoS,2L),4,c(0,1))[2,(Dim:(Dim+1))])
            Res2 <- cbind(SpI,SpS)
            return(list(Res1,Res2))
        }
        degree <- 3
        dbase <- 1
    }

    # Creating the points and weights of Gauss-Legendre quadrature (used when base is "exp.bs" and degree in c(2:3) or when base is "exp.ns")
    gl <- gauss.quad(n=n.gleg,kind="legendre")
    gln <- gl$nodes
    lglw <- log(gl$weights)

    # Formatting the data
    time.new.0 <- NULL
    fix.new <- model.matrix(FormulaF,data=Test)
    names.fix <- colnames(fix.new)[-1]
    nph.new <- model.matrix(FormulaN,data=Test)
    nbtd <- dim(nph.new)[2]
    names.nph <- colnames(nph.new)[-1]

    # Creation of the different objects necessary to compute the hazards

    # Weibull hazard
    if (base=="weibull"){
        time.cat.0 <- NULL
        time.cat <- NULL
        vec.knots <- NULL
        int.knots <- NULL
        MatK <- NULL
        NsAdj <- list(NULL,NULL)
        time.new <- log(time.pts)
        n.td.base <- 1
        n.ntd <- dim(fix.new)[2]
        if (n.ntd>1){
            which.ntd <- c(1,2+1:(n.ntd-1))
        }
        else {
            which.ntd <- 1
        }
        n.td.nph <- length(names.nph)
        if (n.td.nph>0){
            which.td <- c(2,(2+(n.ntd-1))+1:n.td.nph)
        }
        else {
            which.td <- 2
        }
        if (!is.null(cluster)){
            fix.new <- as.data.frame(fix.new)
            fix.new$cluster <- 1
            fix.new <- as.matrix(fix.new)
            which.ntd <- c(which.ntd,n.par)
        }
        fix.new <- fix.new[,-1,drop=FALSE]
        nph.new <- nph.new[,-1,drop=FALSE]
    }

    # Hazard modelled by the exponential of a B-spline
    else if (base%in%c("exp.bs","exp.ns","pw.cst")){

        cuts <- c(0,knots,max.time)
        time.cat.lab <- cut(time.pts,breaks=cuts,include.lowest=TRUE)
        time.cat <- as.numeric(time.cat.lab)-1
        time.cat.0 <- NULL

        if (base=="pw.cst"){
            time.new <- time.pts-cuts[time.cat+1]
            vec.knots <- NULL
            int.knots <- NULL
            MatK <- c(knots,BoS)-c(0,knots)
            NsAdj <- list(NULL,NULL)
            n.td.base <- length(knots)+1

            # For non time-dependent effects
            fix.new <- fix.new[,-which(colnames(fix.new)%in%colnames(nph.new)),drop=FALSE]
            names.fix <- colnames(fix.new)
            if (!is.null(names.fix)){
                n.ntd <- dim(fix.new)[2]
            }
            else {
                n.ntd <- 0
                fix.new <- 0
            }
            if (n.ntd>0){
                which.ntd <- c(n.td.base+1:n.ntd)
            }
            else {
                which.ntd <- NULL
            }
            n.inter <- 0
        }
        else if (base%in%c("exp.bs","exp.ns")){
            time.new <- time.pts
            vec.knots <- c(rep(BoI,degree),knots,rep(BoS,degree))
            ltk <- length(vec.knots)
            MatK <- TransfDeg(vec.knots,degree)
            NsAdj <- NsAdjust(c(BoI,vec.knots,BoS),BoI,BoS)
            if (base=="exp.bs") int.knots <- cuts
            else {
                BOI <- NULL
                BOS <- NULL
                firstK <- 0
                if (BoI>0){
                    BOI <- BoI
                    MatK <- cbind(0,MatK)
                    vec.knots <- c(BoI,vec.knots)
                    firstK <- 1
                }
                if (BoS<max.time){
                    BOS <- BoS
                    MatK <- cbind(MatK,0)
                    vec.knots <- c(vec.knots,BoS)
                }
                int.knots <- c(0,BOI,knots,BOS,max.time)
                time.cat.lab <- cut(time.pts,breaks=int.knots,include.lowest=TRUE)
                time.cat <- as.numeric(time.cat.lab)-1
                degree <- c(degree,ltk,firstK) # ltk and firstK are passed to the HazardInt function through the 'degree' argument which is not used with NS
            }
            n.td.base <- dbase + length(knots)

            # For non time-dependent effects
            n.ntd <- dim(fix.new)[2]
            if (n.ntd>1){
                which.ntd <- c(1,(n.td.base+1)+1:(n.ntd-1))
            }
            else {
                which.ntd <- 1
            }
            n.inter <- 1
        }

        # For time-dependent effects
        n.td.nph <- length(names.nph)*n.td.base
        if (n.td.nph>0){
            which.td <- c(n.inter+1:n.td.base,(n.td.base+n.ntd)+1:n.td.nph)
        }
        else {
            which.td <- c(n.inter+1:n.td.base)
        }

        if (!is.null(cluster)){
            fix.new <- as.data.frame(fix.new)
            fix.new$cluster <- 1
            fix.new <- as.matrix(fix.new)
            which.ntd <- c(which.ntd,n.par)
        }

        fix.new <- t(fix.new) # Matrix of prop effects has to be transposed for use by the Int/Delta functions
        nph.new <- t(nph.new) # Matrix of time-dep effects has to be transposed for use by the Int/Delta functions
    }

    nb.par <- length(c(which.ntd,which.td))
    temp.H <- HazardInt(x0=time.new.0,x=time.new,nph=nph.new,timecat0=time.cat.0,timecat=time.cat,fixobs=fix.new,param=coef[which.td],paramf=coef[which.ntd],deg=degree,n=gln,lw=lglw,matk=MatK,totk=vec.knots,intk=int.knots,nsadj1=NsAdj[[1]],nsadj2=NsAdj[[2]])
    lambda <- exp(temp.H$LogHaz)
    Surv <- exp(-exp(temp.H$LogCum))

    # Variance and CI estimation
    if (conf.int=="delta"){
        varcov <- vcov[c(which.ntd,which.td),c(which.ntd,which.td),drop=FALSE]
        temp.V <- DeltaFct(x0=time.new.0,x=time.new,nph=nph.new,timecat0=time.cat.0,timecat=time.cat,fixobs=fix.new,paramt=c(coef[which.ntd],coef[which.td]),deg=degree,n=gln,lw=lglw,matk=MatK,totk=vec.knots,intk=int.knots,nsadj1=NsAdj[[1]],nsadj2=NsAdj[[2]],varcov=varcov,as.numeric(include.gradient))
        Var.Log.Haz <- temp.V$VarLogHaz
        Var.Log.Cum <- temp.V$VarLogCum
        varcov <- varcov[order(c(which.ntd,which.td)),order(c(which.ntd,which.td)),drop=FALSE]
        if (include.gradient==TRUE){
            Grad.LH <- matrix(temp.V$GradLogHaz,nb.time.pts,nb.par)[,order(c(which.ntd,which.td)),drop=FALSE]
            Grad.LC <- matrix(temp.V$GradLogCum,nb.time.pts,nb.par)[,order(c(which.ntd,which.td)),drop=FALSE]
            colnames(Grad.LH) <- colnames(varcov)
            colnames(Grad.LC) <- colnames(varcov)
        }
        alpha <- (1-level)/2
        if (delta.type.h=="log"){
            BInf1 <- exp(temp.H$LogHaz + qnorm(alpha)*sqrt(Var.Log.Haz))
            BSup1 <- exp(temp.H$LogHaz + qnorm(1-alpha)*sqrt(Var.Log.Haz))
        }
        if (delta.type.h=="plain"){
            Var.Haz <- exp(2*temp.H$LogHaz)*Var.Log.Haz
            BInf1 <- exp(temp.H$LogHaz) + qnorm(alpha)*sqrt(Var.Haz)
            BSup1 <- exp(temp.H$LogHaz) + qnorm(1-alpha)*sqrt(Var.Haz)
        }
        if (delta.type.s=="log-log"){
            BSup2 <- exp(-exp(temp.H$LogCum + qnorm(alpha)*sqrt(Var.Log.Cum)))
            BInf2 <- exp(-exp(temp.H$LogCum + qnorm(1-alpha)*sqrt(Var.Log.Cum)))
        }
        if (delta.type.s=="log"){
            Var.Cum <- exp(2*temp.H$LogCum)*Var.Log.Cum
            BSup2 <- exp(-exp(temp.H$LogCum) + qnorm(alpha)*sqrt(Var.Cum))
            BInf2 <- exp(-exp(temp.H$LogCum) + qnorm(1-alpha)*sqrt(Var.Cum))
        }
        if (delta.type.s=="plain"){
            Var.Surv <- exp(2*(temp.H$LogCum-exp(temp.H$LogCum)))*Var.Log.Cum
            BInf2 <- exp(-exp(temp.H$LogCum)) + qnorm(alpha)*sqrt(Var.Surv)
            BSup2 <- exp(-exp(temp.H$LogCum)) + qnorm(1-alpha)*sqrt(Var.Surv)
        }
    }
    else if (conf.int=="simul"){
        varcov <- vcov
        Res1 <- matrix(0,nb.time.pts,nb.sim)
        Res2 <- matrix(0,nb.time.pts,nb.sim)
        Coef <- mvrnorm(nb.sim,mu=coef,Sigma=varcov)
        for (i in 1:nb.sim){
            p.td <- Coef[i,which.td]
            p.ntd <- Coef[i,which.ntd]
            temp.H <- HazardInt(x0=time.new.0,x=time.new,nph=nph.new,timecat0=time.cat.0,timecat=time.cat,fixobs=fix.new,param=p.td,paramf=p.ntd,deg=degree,n=gln,lw=lglw,matk=MatK,totk=vec.knots,intk=int.knots,nsadj1=NsAdj[[1]],nsadj2=NsAdj[[2]])
            Res1[,i] <- exp(temp.H$LogHaz)
            Res2[,i] <- exp(-exp(temp.H$LogCum))
        }
        Var.Log.Haz <- apply(Res1,1,FUN=var)
        Var.Log.Cum <- apply(Res2,1,FUN=var)
        BInf1 <- apply(Res1,1,quantile,prob=0.025)
        BSup1 <- apply(Res1,1,quantile,prob=0.975)
        BInf2 <- apply(Res2,1,quantile,prob=0.025)
        BSup2 <- apply(Res2,1,quantile,prob=0.975)
    }

    variances <- NA
    df.ct <- NA
    if (conf.int!="none"){
        res1 <- data.frame(hazard=lambda,hazard.inf=BInf1,hazard.sup=BSup1,
                           surv=Surv,surv.inf=BInf2,surv.sup=BSup2)
        if (conf.int=="delta"){
            variances <- data.frame(var.log.haz=Var.Log.Haz,var.log.cum=Var.Log.Cum)
            df.ct <- data.frame(hazard=delta.type.h,surv=delta.type.s)
        }
    }
    else {
        res1 <- data.frame(hazard=lambda,surv=Surv)
    }

    if (!is.null(cluster)){
        Test <- cbind(Test,cluster)
    }
    res.PS <- list(call=call,
                   results=cbind(time.pts,Test,res1),
                   variances=variances,
                   grad.loghaz=Grad.LH,
                   grad.logcum=Grad.LC,
                   vcov=varcov,
                   type=typepred,
                   ci.method=conf.int,
                   level=level,
                   delta.type=df.ct,
                   nb.sim=ifelse(conf.int=="simul",nb.sim,NA))
    res.PS <- res.PS[!is.na(res.PS)]

    class(res.PS) <- "predMexhaz"
    res.PS
}

mexhazEgh <- function(formula,data,expected=NULL,base="weibull",degree=3,knots=NULL,bound=NULL,n.gleg=20,init=NULL,random=NULL,n.aghq=10,verbose=0,iterlim=10000,print.level=1,gradtol=1e-8,mode="fit",...){

    time0 <- as.numeric(proc.time()[3])
    FALCenv <- environment()

    call <- match.call()
    mf <- match.call(expand.dots=FALSE)
    m <- match(c("formula","data"),names(mf),0L)
    if (m[1]==0){
        stop("The 'formula' argument is required...")
    }
    if (m[2]==0){
        stop("The 'data' argument is required...")
    }
    name.data <- paste(substitute(data),sep="")
    if (length(name.data)>1){
        name.data <- name.data[2]
    }
    if (base == "exp.bs" & !degree %in% c(1:3)) {
        stop("This function can only be used to estimate log-hazards described by B-splines of degree 1 to 3...")
    }
    if (!is.null(random) & !is.null(expected)){
        stop("mexhazEgh cannot be used for estimating an excess hazard model with a frailty. Please use the mexhazStd function instead...")
    }

    ## Function that controls what is printed during the optimisation procedure
    if (verbose>0){
        verbose.ll <- function(){
            if (!((iv <- parent.neval/verbose)-floor(iv))){
                time1 <- as.numeric(proc.time()[3])
                print(data.frame(Eval=parent.neval,LogLik=-parent.logLik,Time=time1-time0,row.names=""))
                cat("Param\n")
                print(round(parent.param,4))
                cat("\n")
            }
        }
    }
    else {
        verbose.ll <- function(){}
    }

    ## Functions for computing the part of the B-spline bases that depends only on the knots and can therefore be calculated only once at the beginning of the function (used in combination with the IntHazard function to estimate the hazard and the cumulative hazard)

    ## For B-splines (also used for Natural Cubic Splines)
    TransfDeg <- function(vec.knots,deg){
        if (deg==1){
            Le <- length(vec.knots)
            Res <- 1/(vec.knots[2:Le]-vec.knots[1:(Le-1)])
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

    ## For Natural Cubic Splines
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

    ## Creating the points and weights of Gauss-Legendre quadrature (used when base is "exp.bs" and degree in c(2:3) or when base is "exp.ns")
    if (n.gleg<=0 | round(n.gleg,0)!=n.gleg){
        stop("The 'n.gleg' argument must be a strictly positive integer.")
    }
    gl <- gauss.quad(n=n.gleg,kind="legendre")
    gln <- gl$nodes
    lglw <- log(gl$weights)

    MyProd <- function(mat1,mat2){
        lem1 <- dim(mat1)[2]
        lem2 <- dim(mat2)[2]
        matrix(unlist(sapply(1:lem1,function(i){mat1[,i]*mat2[,i:lem2,drop=FALSE]},simplify=FALSE)),dim(mat1)[1],lem1*(lem2-0.5*(lem1-1)))
    }

    Reorder <- function(vec,mat,idx){
        mat[idx] <- vec
        mat <- t(mat)
        mat[idx] <- vec
        return(mat)
    }

    ## Data preparation
    mf <- mf[c(1L,m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- "na.pass"
    mf[[1L]] <- quote(stats::model.frame)
    tot.formula <- terms(formula,data=data,specials="nph")
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
    mf$formula <- update(tot.formula,paste(".~.",nTerm2,sep="-"))
    FormulaF <- update(tot.formula,paste("NULL~.",nTerm2,sep="-"))
    Xlevel.formula <- terms(FormulaF,data=data)
    m <- eval(mf,parent.frame())
    if (nrow(m)==0){
        stop("No non-missing observations...")
    }

    Y <- model.extract(m,"response")
    if (!inherits(Y,"Surv")){
        stop("Response must be a Surv() object...")
    }
    Survtype <- attr(Y, "type")
    if (ncol(Y)==2){
        if (Survtype!="right"){
            stop(paste("mexhaz does not support \"", Survtype, "\" type of censoring with (0, time] survival data",
                       sep = ""))
        }
        time.obs <- Y[,1,drop=FALSE]   # Follow-up time
        status.obs <- Y[,2]   # Status variable
    }
    else if (ncol(Y)==3){
        if (Survtype!="counting"){
            stop(paste("mexhaz does not support \"", Survtype, "\" type of censoring with (time, time2] survival data",
                       sep = ""))
        }
        time.obs <- Y[,1:2]   # Entry time / Follow-up time
        status.obs <- Y[,3]   # Status variable
    }

    if (Survtype!="counting"){
        if (!is.null(expected)){
            if (base=="weibull"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HGHAggr_WbRx,x,nph,fixobs,statobs,lambdaobs,nbyclust,param,paramf)
                }
            }
            else if (base=="pw.cst"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HGHAggr_PwRx,x,nph,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,as.double(matk))
                }
            }
            else if (base=="exp.bs"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HGHAggr_BsRx,x,nph,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
                }
            }
            else if (base=="exp.ns"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HGHAggr_NsRx,x,nph,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2))
                }
            }
        }
        else {
            if (base=="weibull"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call("HGHAggr_WbR",x,nph,fixobs,statobs,nbyclust,param,paramf)
                }
            }
            else if (base=="pw.cst"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call("HGHAggr_PwR",x,nph,timecat,fixobs,statobs,nbyclust,param,paramf,as.double(matk))
                }
            }
            else if (base=="exp.bs"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call("HGHAggr_BsR",x,nph,timecat,fixobs,statobs,nbyclust,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
                }
            }
            else if (base=="exp.ns"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call("HGHAggr_NsR",x,nph,timecat,fixobs,statobs,nbyclust,param,paramf,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2))
                }
            }
        }
        HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
            temp <- HazGradHess(x0,x,t(nph),timecat0,timecat,t(fixobs),statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2)
            nclust <- length(nbyclust)
            lp <- length(paramf)+length(param)
            Result <- list(LogHaz=temp$LogHaz, HazCum=temp$HazCum, Grad.loglambda=matrix(temp$GradLogHaz,nclust,lp),
                           Grad.Lambda=matrix(temp$GradCum,nclust,lp), Hess.loglambda=matrix(temp$HessLHaz,nclust,0.5*lp*(lp+1)),
                           Hess.Lambda=matrix(temp$HessCum,nclust,0.5*lp*(lp+1)), Test=temp$Test)
            return(Result)
        }
    }
    else {
        if (!is.null(expected)){
            if (base=="weibull"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HGHAggr_WbLx,x0,x,nph,fixobs,statobs,lambdaobs,nbyclust,param,paramf)
                }
            }
            else if (base=="pw.cst"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HGHAggr_PwLx,x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,as.double(matk))
                }
            }
            else if (base=="exp.bs"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HGHAggr_BsLx,x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
                }
            }
            else if (base=="exp.ns"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HGHAggr_NsLx,x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2))
                }
            }
        }
        else {
            if (base=="weibull"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call("HGHAggr_WbL",x0,x,nph,fixobs,statobs,nbyclust,param,paramf)
                }
            }
            else if (base=="pw.cst"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call("HGHAggr_PwL",x0,x,nph,timecat0,timecat,fixobs,statobs,nbyclust,param,paramf,as.double(matk))
                }
            }
            else if (base=="exp.bs"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call("HGHAggr_BsL",x0,x,nph,timecat0,timecat,fixobs,statobs,nbyclust,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
                }
            }
            else if (base=="exp.ns"){
                HazGradHess <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call("HGHAggr_NsL",x0,x,nph,timecat0,timecat,fixobs,statobs,nbyclust,param,paramf,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2))
                }
            }
        }
        HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
            temp <- HazGradHess(x0,x,t(nph),timecat0,timecat,t(fixobs),statobs,lambdaobs,nbyclust,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2)
            nclust <- length(nbyclust)
            lp <- length(paramf)+length(param)
            Result <- list(LogHaz=temp$LogHaz, HazCum0=temp$HazCum0, HazCum=temp$HazCum, Grad.loglambda=matrix(temp$GradLogHaz,nclust,lp),
                           Grad.Lambda=matrix(temp$GradCum,nclust,lp), Grad.Lambda0=matrix(temp$GradCum0,nclust,lp),
                           Hess.loglambda=matrix(temp$HessLHaz,nclust,0.5*lp*(lp+1)), Hess.Lambda=matrix(temp$HessCum,nclust,0.5*lp*(lp+1)),
                           Hess.Lambda0=matrix(temp$HessCum0,nclust,0.5*lp*(lp+1)), Test=temp$Test)
            return(Result)
        }
    }

    ApproxLamb <- function(x){
        l2 <- log(x)
        x-l2+l2/x+(l2*(l2-2)/(2*x^2))+(l2*(2*l2^2-9*l2+6)/(6*x^3))+(l2*(3*l2^3-22*l2^2+36*l2-12))/(12*x^4)+(l2*(12*l2^4-125*l2^3+350*l2^2-300*l2+60))/(60*x^5)
    }

    data.fix <- model.frame(FormulaF,data=data,na.action=na.pass)
    data.nph <- model.frame(FormulaN,data=data,na.action=na.pass)
    n.obs.tot <- dim(data)[1] # Total number of observations

    if (!is.null(expected)) {
        lambda.pop <- data[, expected]
        if (min(lambda.pop, na.rm = TRUE) < 0) {
            stop("The expected hazard for some observations is negative!")
        }
        withExp <- 1
    }
    else {
        lambda.pop <- rep(0, n.obs.tot)
        withExp <- 0
    }

    ## Remove observations containing missing values and perform several formatting operations on the data
    if (!is.null(random)){
        random.obs <- data[,random]
        Idx.NA.Rdm <- which(is.na(random.obs))
        if (length(Idx.NA.Rdm)>0){
            warning("Cluster information was missing for some observations. These observations were consequently removed...")
            random.obs <- random.obs[-Idx.NA.Rdm]
            time.obs <- time.obs[-Idx.NA.Rdm,,drop=FALSE]
            status.obs <- status.obs[-Idx.NA.Rdm]
            data.fix <- data.fix[-Idx.NA.Rdm,,drop=FALSE]
            data.nph <- data.nph[-Idx.NA.Rdm,,drop=FALSE]
        }
        random.obs <- as.factor(random.obs)
        clust <- levels(random.obs)
        IdxRnd <- order(random.obs)
        ## The dataset MUST be ordered by the levels of the clustering variable
        random.obs <- random.obs[IdxRnd]
        time.obs <- time.obs[IdxRnd,,drop=FALSE]
        status.obs <- status.obs[IdxRnd]
        lambda.pop <- lambda.pop[IdxRnd]
        data.fix <- data.fix[IdxRnd,,drop=FALSE]
        data.nph <- data.nph[IdxRnd,,drop=FALSE]
    }
    TotData <- cbind(time.obs,status.obs,lambda.pop,data.fix,data.nph)
    Xlevels <- .getXlevels(Xlevel.formula,TotData)
    Idx.Non.NA <- which(apply(TotData,1,function(vec){sum(is.na(vec))})==0)
    if (length(Idx.Non.NA)<n.obs.tot){
        warning("Covariables information was missing for some observations. These observations were consequently removed...")
    }
    time.obs <- time.obs[Idx.Non.NA,,drop=FALSE]
    n.obs <- dim(time.obs)[1]
    if (n.obs==0){
        stop("No non-missing values for some covariables...")
    }
    if (Survtype=="right"){
        time.obs.0 <- NULL
        Idx.Time.Neg <- which(time.obs<0)
        if (length(Idx.Time.Neg)>0){
            stop("Some observations have a negative follow-up time...")
        }
        Idx.Time.0 <- which(time.obs==0)
        if (length(Idx.Time.0)>0){
            warning("Some observations had a follow-up time of length 0. A value of 1/730.5 (half a day) was substituted. Please check to see if it is appropriate or deal with 0 follow-up time values before using the mexhaz function.")
        }
        time.obs[Idx.Time.0] <- 1/730.5
        py.obs <- sum(time.obs)
    }
    if (Survtype=="counting"){
        time.obs.0 <- time.obs[,1]
        time.obs <- time.obs[,2]
        Idx.Time.Neg1 <- which(time.obs.0<0)
        if (length(Idx.Time.Neg1)>0){
            stop("Some observations have a negative entry time...")
        }
        Idx.Time.Neg2 <- which(time.obs<0)
        if (length(Idx.Time.Neg2)>0){
            stop("Some observations have a negative exit time...")
        }
        py.obs <- sum(time.obs-time.obs.0)
    }
    max.time <- max(time.obs)
        TestBo <- 0
    if (!is.null(bound) & base%in%c("exp.bs","exp.ns")){
        if (!is.numeric(bound) | length(bound)!=2 | bound[1]==bound[2]){
            warning("The 'bound' argument should be a vector of two non-equal numeric values. Consequently, the boundary knots were assigned the default values (0,max.time).")
            TestBo <- 1
        }
        else {
            if (bound[2]<bound[1]){
                bound <- bound[order(bound)]
                warning("The 'bound' argument should be a vector of two ordered numeric values. Consequently, the values given were re-ordered.")
            }
            if (base=="exp.bs" & (bound[1]>0 | bound[2]<max.time)){
                warning("When used with the 'exp.bs' hazard, bound[1] should not be greater than 0 and bound[2] should not be lower than max.time. Consequently, the boundary knots were assigned the default values (0,max.time).")
                TestBo <- 1
            }
        }
    }
    else {
        TestBo <- 1
    }
    if (TestBo==1){
        BoI <- 0
        BoS <- max.time
    }
    else {
        BoI <- bound[1]
        BoS <- bound[2]
    }

    status.obs <- status.obs[Idx.Non.NA]
    n.events <- sum(status.obs)
    data.fix <- data.fix[Idx.Non.NA,,drop=FALSE]
    fix.obs <- model.matrix(FormulaF,data=data.fix,drop.unused.levels=TRUE)
    names.fix <- colnames(fix.obs)[-1]
    data.nph <- data.nph[Idx.Non.NA,,drop=FALSE]
    nph.obs <- model.matrix(FormulaN,data=data.nph,drop.unused.levels=TRUE)
    nbtd <- dim(nph.obs)[2]
    names.nph <- colnames(nph.obs)[-1]
    if (sum(!names.nph%in%names.fix)){
        stop("Some variables in the 'nph()' sub-formula are not included in the formula...")
    }
    lambda.pop <- lambda.pop[Idx.Non.NA]

    ## Remove unnecessary objects
    rm(TotData,data.fix,data.nph)

    ## Creation of the different objects necessary to compute the hazards

    ## Weibull hazard
    if (base=="weibull"){
        time.cat.0 <- NULL
        time.cat <- NULL
        degree <- NA
        vec.knots <- NULL
        int.knots <- NULL
        MatK <- NULL
        NsAdj <- list(NULL,NULL)
        n.td.base <- 1
        n.ntd <- dim(fix.obs)[2]
        if (n.ntd>1){
            which.ntd <- c(1,2+1:(n.ntd-1))
        }
        else {
            which.ntd <- 1
        }
        n.td.nph <- length(names.nph)
        if (n.td.nph>0){
            names.nph <- paste("Rho",names.nph,sep="*")
            which.td <- c(2,(2+(n.ntd-1))+1:n.td.nph)
        }
        else {
            which.td <- 2
        }
        param.names <- c("logLambda","logRho",names.fix,names.nph)
        n.par.fix <- n.td.base+n.ntd+n.td.nph
        param.init <- rep(0,n.par.fix)
    }

    ## Hazard modelled by the exponential of a B-spline / restricted cubic B-spline / Piecewise constant
    else if (base%in%c("exp.bs","exp.ns","pw.cst")){
        ## Baseline hazard-related objects
        if (!is.null(knots)){
            if (min(knots,na.rm=TRUE)<=0 | max(knots,na.rm=TRUE)>=max.time | is.na(sum(knots))){
                stop("The 'knots' argument should be a vector of values strictly between 0 and max.time, the maximum follow-up time observed in the dataset")
            }
            else {
                if (sum(abs(knots-knots[order(knots)]))>0){
                    warning("The knots were not given in increasing order. They were consequently re-ordered.")
                    knots <- knots[order(knots)]
                }
                if (length(unique(knots))<length(knots)){
                    warning("There were duplicate values in the vector of knots. Duplicate values were removed.")
                    knots <- unique(knots)
                }
            }
        }
        cuts <- c(0,knots,max.time)
        time.cat.lab <- cut(time.obs,breaks=cuts,include.lowest=TRUE)
        time.cat <- as.numeric(time.cat.lab)-1
        if (base=="pw.cst"){
            degree <- 0
            time.obs <- time.obs-cuts[time.cat+1]
            if (Survtype=="counting"){
                time.cat.0 <- as.numeric(cut(time.obs.0,breaks=cuts,include.lowest=TRUE))-1
                time.obs.0 <- time.obs.0-cuts[time.cat.0+1]
            }
            else if (Survtype=="right"){
                time.cat.0 <- NULL
            }
            vec.knots <- NULL
            int.knots <- NULL
            MatK <- c(knots,BoS)-c(BoI,knots)
            NsAdj <- list(NULL,NULL)
            n.td.base <- length(knots)+1
            names.base <- levels(time.cat.lab)

            ## For non time-dependent effects
            fix.obs <- fix.obs[,-which(colnames(fix.obs)%in%colnames(nph.obs)),drop=FALSE]
            names.fix <- colnames(fix.obs)
            if (!is.null(names.fix)){
                n.ntd <- dim(fix.obs)[2]
            }
            else {
                n.ntd <- 0
                fix.obs <- 0
            }
            if (n.ntd>0){
                which.ntd <- c(n.td.base+1:n.ntd)
            }
            else {
                which.ntd <- NULL
            }
            intercept <- NULL
            n.inter <- 0
        }
        else if (base%in%c("exp.bs","exp.ns")){
            if (Survtype=="counting"){
                time.cat.0 <- as.numeric(cut(time.obs.0,breaks=cuts,include.lowest=TRUE))-1
            }
            else if (Survtype=="right"){
                time.cat.0 <- NULL
            }
            vec.knots <- c(rep(BoI,degree),knots,rep(BoS,degree))
            ltk <- length(vec.knots)
            MatK <- TransfDeg(vec.knots,degree)
            NsAdj <- NsAdjust(c(BoI,vec.knots,BoS),BoI,BoS)
            n.td.base <- dbase + length(knots)
            names.base <- paste(ifelse(base=="exp.bs","BS","NS"),degree,".",1:n.td.base,sep = "")

            if (base=="exp.bs") {
                int.knots <- cuts
            }
            else {
                BOI <- NULL
                BOS <- NULL
                firstK <- 1
                if (BoI>0){
                    BOI <- BoI
                    MatK <- cbind(0,MatK)
                    vec.knots <- c(BoI,vec.knots)
                    firstK <- 0
                }
                if (BoS<max.time){
                    BOS <- BoS
                    MatK <- cbind(MatK,0)
                    vec.knots <- c(vec.knots,BoS)
                }
                int.knots <- c(0,BOI,knots,BOS,max.time)
                time.cat.lab <- cut(time.obs,breaks=int.knots,include.lowest=TRUE)
                time.cat <- as.numeric(time.cat.lab)-1
                degree <- c(degree,ltk,firstK) # ltk and firstK are passed to the HazardInt function through the 'degree' argument which is not used with NS
            }

            ## For non time-dependent effects
            n.ntd <- dim(fix.obs)[2]
            if (n.ntd>1){
                which.ntd <- c(1,(n.td.base+1)+1:(n.ntd-1))
            }
            else {
                which.ntd <- 1
            }
            intercept <- "Intercept"
            n.inter <- 1
        }

        ## For time-dependent effects
        n.td.nph <- length(names.nph)*n.td.base
        if (n.td.nph>0){
            which.td <- c(n.inter+1:n.td.base,(n.td.base+n.ntd)+1:n.td.nph)
        }
        else {
            which.td <- c(n.inter+1:n.td.base)
        }
        names.nph <- unlist(sapply(names.nph,function(x){paste(x,names.base,sep="*")}))

        param.names <- c(intercept,names.base,names.fix,names.nph)
        n.par.fix <- n.td.base+n.ntd+n.td.nph
        param.init <- rep(0,n.par.fix)
    }

    ## Creation of objects related to the random effect
    n.clust <- 1
    n.rand <- 0
    n.par <- n.td.base + n.ntd + n.td.nph
    which.rdm <- NULL
    if (!is.null(random)){

        n.rand <- 1
        n.par <- n.par + 1
        which.rdm <- n.par
        random.obs <- random.obs[Idx.Non.NA]
        status.one <- which(status.obs==1)
        n.clust <- length(clust) # Number of clusters
        n.by.clust <- as.vector(table(random.obs)) # Size of each cluster
        nbcd <- as.vector(table(random.obs[status.one]))
        Idx.clust.no.obs <- which(n.by.clust==0)
        listIdx <- sapply(levels(random.obs),function(x){which(random.obs==x)},simplify=FALSE)
        if (Survtype=="counting"){
            listIdx0 <- sapply(levels(random.obs),function(x){which(random.obs==x & time.obs.0!=0)},simplify=FALSE)
        }
        if (length(Idx.clust.no.obs)>0){
            warning("Some clusters had no valid observations and were consequently removed...")
            cat("The following clusters were removed:\n")
            print(names(table(random.obs))[Idx.clust.no.obs])
            n.by.clust <- n.by.clust[-Idx.clust.no.obs]
            nbcd <- nbcd[-Idx.clust.no.obs]
            clust <- clust[-Idx.clust.no.obs]
            n.clust <- length(clust)
        }
        nbcd2 <- nbcd^2
        param.names <- c(param.names,paste(random," [log(sd)]",sep=""))

        ## Creation of the points and weights of Gauss-Hermite quadrature
        if (n.aghq<=0 | round(n.aghq,0)!=n.aghq){
            stop("The 'n.aghq' argument must be a strictly positive integer.")
        }
        gq <- gauss.quad(n=n.aghq,kind="hermite")
        x.H <- gq$nodes
        rho.H <- gq$weights

        SumB <- function(x,v){
            temp <- sqrt(2*v/(1+x))*x.H
            rho.H%*%exp(-x/v*(exp(temp)-0.5*(1+(temp+1)^2)))
        }

        SumHB <- function(x,v){
            temp <- sqrt(2*v/(1+x))*x.H
            g <- exp(temp)-0.5*(1+(temp+1)^2)
            gp <- x.H*(exp(temp)-(temp+1))
            gp2 <- x.H^2*(exp(temp)-1)
            exp.tg <- exp(-x/v*g)
            rho.H%*%cbind(g*exp.tg,gp*exp.tg,gp2*exp.tg,g^2*exp.tg,g*gp*exp.tg,gp^2*exp.tg)
        }

    }
    else {
        n.clust <- 1
        n.by.clust <- n.obs
    }

    ## Initial values
    Order <- order(c(which.ntd,which.td,which.rdm))
    if (!is.null(init)){
        if (length(init)!=n.par){
            print(data.frame(Time.Dep.baseline=n.td.base,Non.TD.Effects=n.ntd,TD.Effects=n.td.nph,Random.Effect=n.rand))
            stop("Wrong number of init parameters",call.=FALSE)
        }
        Init <- 1
    }
    else {
        Init <- 0
        init <- c(param.init,rep(0.1,n.rand))
        init[1] <- 0.5*log(n.events/py.obs)
    }

    ## Creation of variables that will be modified by the Hazard and LL.Tot functions
    parent.neval <- 0
    parent.param <- init
    parent.logLik <- -.Machine$double.xmax

    ## For the Hessian
    Mat <- matrix(0,n.par,n.par)
    IdxM <- which(lower.tri(Mat,diag=TRUE))
    if (!is.null(random)){
        IdxM2 <- n.par*1:n.par
        IdxM1 <- IdxM[!IdxM%in%IdxM2]
        IdxM <- c(IdxM1,IdxM2)
    }

    ## Function that actually computes the log-likelihood
    LL.Tot <- function(p.LT,mu.hat.LT=0,env=FALCenv){

        p.td <- p.LT[which.td]
        p.ntd <- p.LT[which.ntd]
        temp.H <- HazardInt(x0=time.obs.0,x=time.obs,nph=nph.obs,timecat0=time.cat.0,timecat=time.cat,fixobs=fix.obs,statobs=status.obs,lambdaobs=lambda.pop,nbyclust=n.by.clust,param=p.td,paramf=p.ntd,deg=degree,n=gln,lw=lglw,matk=MatK,totk=vec.knots,intk=int.knots,nsadj1=NsAdj[[1]],nsadj2=NsAdj[[2]])
        if (temp.H$Test){
            res.LT <- .Machine$double.xmax
        }
        else {
            if (!is.null(random)){

                lvar.w <- 2*p.LT[which.rdm]
                var.w <- exp(lvar.w)
                sd.w <- exp(p.LT[which.rdm])

                ## Log-Likelihood

                Hi <- temp.H$HazCum
                lhi <- temp.H$LogHaz
                lwNu <- log(Hi)+lvar.w+var.w*nbcd
                Nu <- lambertW0(exp(lwNu))
                IdxNu <- which(Nu==Inf)
                if (length(IdxNu>0)){
                    Nu[IdxNu] <- ApproxLamb(lwNu[IdxNu])
                }
                SnB <- sapply(Nu,SumB,v=var.w)
                LogLik <- -sum(lhi-0.5*(log(pi*(1+Nu))+Nu*(Nu+2)/var.w-nbcd2*var.w)+log(SnB))

                ## Gradient

                dlhi.db <- temp.H$Grad.loglambda
                dHi.db <- temp.H$Grad.Lambda

                dNu.db <- dHi.db*Nu/((1+Nu)*Hi)
                dNu.ds <- 2*Nu*(1+nbcd*var.w)/(1+Nu)

                Psi <- sqrt(2*var.w/(1+Nu))
                dPsi.db <- -Psi*dNu.db/(2*(1+Nu))
                dPsi.ds <- -Psi*(dNu.ds/(2*(1+Nu))-1)
                dPsi.di <- cbind(dPsi.db,dPsi.ds,deparse.level=0)

                Theta <- -Nu/var.w
                dTheta.di <- -(1/var.w)*cbind(dNu.db,dNu.ds-2*Nu,deparse.level=0)

                C1 <- (1/(1+Nu)+2*(Nu+1)/var.w)
                dF.db <- dlhi.db-0.5*C1*dNu.db
                dF.ds <- -0.5*C1*dNu.ds+(Nu*(Nu+2)/var.w+nbcd2*var.w)
                dF.di <- cbind(dF.db,dF.ds,deparse.level=0)
                sumHB <- t(sapply(Nu,SumHB,v=var.w))/SnB
                dlSn.di <- dTheta.di*sumHB[,1]+Theta*dPsi.di*sumHB[,2]
                dLL.di <- apply(dF.di+dlSn.di,2,sum)

                ## Hessian

                d2lhi.db2 <- temp.H$Hess.loglambda
                d2Hi.db2 <- temp.H$Hess.Lambda

                NiN <- Nu/(1+Nu)
                TempProd <- MyProd(dHi.db,dHi.db)
                d2Nu.db2 <- (NiN/Hi)*(-(Nu*(Nu+2)/(Hi*(1+Nu)^2))*TempProd+d2Hi.db2)
                d2Nu.db.ds <- 2*NiN*((1+nbcd*var.w)/(Hi*(1+Nu)^2))*dHi.db
                d2Nu.ds2 <- 4*NiN*(((1+nbcd*var.w)/(1+Nu))^2+nbcd*var.w)
                d2Nu.di2 <- cbind(d2Nu.db2,d2Nu.db.ds,d2Nu.ds2,deparse.level=0)

                PiP <- 0.5*Psi/(1+Nu)
                dNu2.db <- TempProd*(Nu/((1+Nu)*Hi))^2
                d2Psi.db2 <- PiP*(1.5*dNu2.db/(1+Nu)-d2Nu.db2)
                d2Psi.db.ds <- PiP*(1.5*dNu.db*dNu.ds/(1+Nu)-d2Nu.db.ds-dNu.db)
                d2Psi.ds2 <- PiP*(1.5*dNu.ds^2/(1+Nu)-d2Nu.ds2-2*dNu.ds+2*(1+Nu))
                d2Psi.di2 <- cbind(d2Psi.db2,d2Psi.db.ds,d2Psi.ds2,deparse.level=0)

                d2Theta.db2 <- -d2Nu.db2/var.w
                d2Theta.db.ds <- (-d2Nu.db.ds+2*dNu.db)/var.w
                d2Theta.ds2 <- (-d2Nu.ds2+4*dNu.ds-4*Nu)/var.w
                d2Theta.di2 <- cbind(d2Theta.db2,d2Theta.db.ds,d2Theta.ds2,deparse.level=0)

                C2 <- (2/var.w-1/(1+Nu)^2)
                d2F.db2 <- d2lhi.db2-0.5*(C1*d2Nu.db2+C2*dNu2.db)
                d2F.db.ds <- -0.5*(C1*d2Nu.db.ds+C2*dNu.db*dNu.ds-4*(Nu+1)/var.w*dNu.db)
                d2F.ds2 <- -0.5*(C1*d2Nu.ds2+C2*dNu.ds^2-8*(Nu+1)/var.w*dNu.ds+4*(Nu*(Nu+2)/var.w-nbcd2*var.w))
                d2F.di2 <- cbind(d2F.db2,d2F.db.ds,d2F.ds2,deparse.level=0)

                dTP.di.PT.di <- MyProd(dTheta.di,dPsi.di)+MyProd(dPsi.di,dTheta.di)
                dPsiPsi.di <- MyProd(dPsi.di,dPsi.di)
                d2lSn.di2 <- d2Theta.di2*sumHB[,1] +
                    (dTP.di.PT.di+Theta*d2Psi.di2)*sumHB[,2] +
                    Theta*dPsiPsi.di*sumHB[,3] +
                    MyProd(dTheta.di,dTheta.di)*sumHB[,4] +
                    Theta*dTP.di.PT.di*sumHB[,5] +
                    Theta^2*dPsiPsi.di*sumHB[,6] -
                    MyProd(dlSn.di,dlSn.di)
                LL.Temp <- apply(d2F.di2+d2lSn.di2,2,sum)
                hess <- Reorder(LL.Temp,Mat,IdxM)

                if (Survtype=="counting"){

                    ## Log-Likelihood

                    Hi0 <- temp.H$HazCum0
                    lwNu0 <- log(Hi0)+lvar.w
                    Nu0 <- lambertW0(exp(lwNu0))
                    IdxNu0 <- which(Nu0==Inf)
                    if (length(IdxNu0>0)){
                        Nu0[IdxNu0] <- ApproxLamb(lwNu0[IdxNu0])
                    }
                    Sn0B <- sapply(Nu0,SumB,v=var.w)
                    LogLik <- -sum(lhi-0.5*(log((1+Nu)/(1+Nu0))+(Nu*(Nu+2)-Nu0*(Nu0+2))/var.w-nbcd2*var.w)+log(SnB)-log(Sn0B))

                    ## Gradient

                    dHi0.db <- temp.H$Grad.Lambda0
                    dNu0.db <- dHi0.db*Nu0/((1+Nu0)*Hi0)
                    dNu0.ds <- 2*Nu0/(1+Nu0)

                    Psi0 <- sqrt(2*var.w/(1+Nu0))
                    dPsi0.db <- -Psi0*dNu0.db/(2*(1+Nu0))
                    dPsi0.ds <- -Psi0*(dNu0.ds/(2*(1+Nu0))-1)
                    dPsi0.di <- cbind(dPsi0.db,dPsi0.ds,deparse.level=0)

                    Theta0 <- -Nu0/var.w
                    dTheta0.di <- -(1/var.w)*cbind(dNu0.db,dNu0.ds-2*Nu0,deparse.level=0)

                    C10 <- (1/(1+Nu0)+2*(Nu0+1)/var.w)
                    dF0.db <- -0.5*C10*dNu0.db
                    dF0.ds <- -0.5*C10*dNu0.ds+Nu0*(Nu0+2)/var.w
                    dF0.di <- cbind(dF0.db,dF0.ds,deparse.level=0)

                    sumHB0 <- t(sapply(Nu0,SumHB,v=var.w))/Sn0B
                    dlSn0.di <- (dTheta0.di*sumHB0[,1] + Theta0*dPsi0.di*sumHB0[,2])
                    dLL0.di <- apply(dF0.di+dlSn0.di,2,sum)
                    dLL.di <- dLL.di-dLL0.di

                    ## Hessian

                    d2Hi0.db2 <- temp.H$Hess.Lambda0

                    NiN0 <- Nu0/(1+Nu0)
                    d2Nu0.db2 <- (NiN0/Hi0)*(-(Nu0*(Nu0+2)/(Hi0*(1+Nu0)^2))*MyProd(dHi0.db,dHi0.db)+d2Hi0.db2)
                    d2Nu0.db.ds <- (2*NiN0/(Hi0*(1+Nu0)^2))*dHi0.db
                    d2Nu0.ds2 <- 4*NiN0/(1+Nu0)^2
                    d2Nu0.di2 <- cbind(d2Nu0.db2,d2Nu0.db.ds,d2Nu0.ds2,deparse.level=0)

                    PiP0 <- 0.5*Psi0/(1+Nu0)
                    dNu02.db <- MyProd(dNu0.db,dNu0.db)
                    d2Psi0.db2 <- PiP0*(1.5*dNu02.db/(1+Nu0)-d2Nu0.db2)
                    d2Psi0.db.ds <- PiP0*(1.5*dNu0.db*dNu0.ds/(1+Nu0)-d2Nu0.db.ds-dNu0.db)
                    d2Psi0.ds2 <- PiP0*(1.5*dNu0.ds^2/(1+Nu0)-d2Nu0.ds2-2*dNu0.ds+2*(1+Nu0))
                    d2Psi0.di2 <- cbind(d2Psi0.db2,d2Psi0.db.ds,d2Psi0.ds2,deparse.level=0)

                    d2Theta0.db2 <- -d2Nu0.db2/var.w
                    d2Theta0.db.ds <- (-d2Nu0.db.ds+2*dNu0.db)/var.w
                    d2Theta0.ds2 <- (-d2Nu0.ds2+4*dNu0.ds-4*Nu0)/var.w
                    d2Theta0.di2 <- cbind(d2Theta0.db2,d2Theta0.db.ds,d2Theta0.ds2,deparse.level=0)

                    C20 <- (2/var.w-1/(1+Nu0)^2)
                    d2F0.db2 <- -0.5*(C10*d2Nu0.db2+C20*dNu02.db)
                    d2F0.db.ds <- -0.5*(C10*d2Nu0.db.ds+C20*dNu0.db*dNu0.ds-4*(Nu0+1)/var.w*dNu0.db)
                    d2F0.ds2 <- -0.5*(C10*d2Nu0.ds2+C20*dNu0.ds^2-8*(Nu0+1)/var.w*dNu0.ds+4*Nu0*(Nu0+2)/var.w)
                    d2F0.di2 <- cbind(d2F0.db2,d2F0.db.ds,d2F0.ds2,deparse.level=0)

                    dTP0.di.PT0.di <- MyProd(dTheta0.di,dPsi0.di)+MyProd(dPsi0.di,dTheta0.di)
                    dPsiPsi0.di <- MyProd(dPsi0.di,dPsi0.di)
                    d2lSn0.di2 <- d2Theta0.di2*sumHB0[,1] +
                        (dTP0.di.PT0.di+Theta0*d2Psi0.di2)*sumHB0[,2] +
                        Theta0*dPsiPsi0.di*sumHB0[,3] +
                        MyProd(dTheta0.di,dTheta0.di)*sumHB0[,4] +
                        Theta0*dTP0.di.PT0.di*sumHB0[,5] +
                        Theta0^2*dPsiPsi0.di*sumHB0[,6] -
                        MyProd(dlSn0.di,dlSn0.di)
                    LL0.Temp <- apply(d2F0.di2+d2lSn0.di2,2,sum)
                    hess <- Reorder(LL.Temp-LL0.Temp,Mat,IdxM)

                }

                if (mu.hat.LT==1){
                    res.LT <- list()
                    res.LT[[1]] <- LogLik
                    res.LT[[2]] <- -dLL.di[Order] # Gradient
                    res.LT[[3]] <- -hess[Order,Order] # Hessian
                    res.LT[[4]] <- var.w*nbcd-Nu # mu.hat
                    res.LT[[5]] <- sqrt(var.w/(1+Nu)) # sigma.hat
                    res.LT[[6]] <- cbind(-dNu.db,2*var.w*nbcd-dNu.ds,deparse.level=0)[,Order] # Gradient mu.hat
                }
                else if (mu.hat.LT==0){
                    res.LT <- LogLik
                    res.LT[is.nan(res.LT) | abs(res.LT)==Inf] <- .Machine$double.xmax
                    attr(res.LT,"gradient") <- -dLL.di[Order]
                    attr(res.LT,"hessian") <- -hess[Order,Order]
                }
            }
            else {
                lhi <- temp.H$LogHaz
                Hi <- temp.H$HazCum

                ## Gradient
                dlhi.db <- temp.H$Grad.loglambda
                dHi.db <- temp.H$Grad.Lambda

                ## Hessian
                d2lhi.db2 <- temp.H$Hess.loglambda
                d2Hi.db2 <- temp.H$Hess.Lambda

                if (Survtype=="counting"){
                    Hi <- Hi - temp.H$HazCum0
                    dHi.db <- dHi.db - temp.H$Grad.Lambda0
                    d2Hi.db2 <- d2Hi.db2 - temp.H$Hess.Lambda0
                }
                if (mu.hat.LT==1){
                    res.LT <- list()
                    res.LT[[1]] <- Hi - lhi
                    res.LT[[2]] <- (dHi.db-dlhi.db)[Order] # Gradient
                    res.LT[[3]] <- Reorder(d2Hi.db2-d2lhi.db2,Mat,IdxM)[Order,Order] # Hessian
                }
                else if (mu.hat.LT==0){
                    res.LT <- Hi - lhi
                    res.LT[is.nan(res.LT) | abs(res.LT)==Inf] <- .Machine$double.xmax
                    attr(res.LT,"gradient") <- (dHi.db-dlhi.db)[Order]
                    attr(res.LT,"hessian") <- Reorder(d2Hi.db2-d2lhi.db2,Mat,IdxM)[Order,Order]
                }
            }
            env$parent.neval <- parent.neval + 1
            if (mu.hat.LT==0){
                env$parent.param <- p.LT + 0 # Necessary...
                env$parent.logLik <- res.LT
                verbose.ll()
            }
       }
        return(res.LT)
    }

    ## Launching the optimisation procedure
    time0 <- as.numeric(proc.time()[3])
    if (mode=="fit"){
        mod.lik <- try(nlm(LL.Tot,init,iterlim=iterlim,print.level=print.level,gradtol=gradtol,check.analyticals=FALSE,...),silent=TRUE)
        if (class(mod.lik)[1]!="try-error"){
            param.fin <- mod.lik$estimate
            code.fin <- mod.lik$code
            loglik <- -mod.lik$minimum
            iterations <- mod.lik$iterations
            names(param.fin) <- param.names
            eval.fin <- LL.Tot(param.fin,mu.hat.LT=1)
            vcov <- try(solve(eval.fin[[3]]),silent=TRUE)
        }
        Flag <- (class(mod.lik)[1]=="try-error")
        if (Flag==0){
            Flag <- (class(vcov)[1]=="try-error")
            if (Flag==0){
                Flag <- (mod.lik$code!=1 | sum(is.na(diag(vcov))| diag(vcov)<0))
            }
        }
        if (Init==0 & Flag){
            Iter <- 1
            while (Iter<9 & Flag){
                Iter <- Iter + 1
                init[1] <- init[1]/Iter
                mod.lik <- try(nlm(LL.Tot,init,iterlim=iterlim,print.level=print.level,gradtol=gradtol,check.analyticals=FALSE,...),silent=TRUE)
                if (class(mod.lik)[1]!="try-error"){
                    param.fin <- mod.lik$estimate
                    code.fin <- mod.lik$code
                    loglik <- -mod.lik$minimum
                    iterations <- mod.lik$iterations
                    names(param.fin) <- param.names
                    eval.fin <- LL.Tot(param.fin,mu.hat.LT=1)
                    vcov <- try(solve(eval.fin[[3]]),silent=TRUE)
                }
                Flag <- (class(mod.lik)[1]=="try-error")
                if (Flag==0){
                    Flag <- (class(vcov)[1]=="try-error")
                    if (Flag==0){
                        Flag <- (mod.lik$code!=1 | sum(is.na(diag(vcov))| diag(vcov)<0))
                    }
                }
            }
        }
        if (Flag){
            Iter <- 1
            while (Iter<9 & Flag){
                init[n.par] <- seq(1,-1,le=9)[Iter]
                mod.lik <- try(nlm(LL.Tot,init,iterlim=iterlim,print.level=print.level,gradtol=gradtol,check.analyticals=FALSE,...),silent=TRUE)
                if (class(mod.lik)[1]!="try-error"){
                    param.fin <- mod.lik$estimate
                    code.fin <- mod.lik$code
                    loglik <- -mod.lik$minimum
                    iterations <- mod.lik$iterations
                    names(param.fin) <- param.names
                    eval.fin <- LL.Tot(param.fin,mu.hat.LT=1)
                    vcov <- try(solve(eval.fin[[3]]),silent=TRUE)
                }
                Flag <- (class(mod.lik)[1]=="try-error")
                if (Flag==0){
                    Flag <- (class(vcov)[1]=="try-error")
                    if (Flag==0){
                        Flag <- (mod.lik$code!=1 | sum(is.na(diag(vcov))| diag(vcov)<0))
                    }
                }
                Iter <- Iter + 1
            }
        }
        if (class(mod.lik)[1]!="try-error"){
            gradient <- eval.fin[[2]]
            hessian <- eval.fin[[3]]
            if (class(vcov)[1]=="try-error"){
                vcov <- matrix(NA,n.par,n.par)
                warning("Unable to invert the Hessian matrix...")
            }
            colnames(vcov) <- rownames(vcov) <- param.names

            ## Empirical Bayes Estimates
            if (!is.null(random)) {
                mu.hat.fin <- eval.fin[[4]]
                sigma.hat.fin <- eval.fin[[5]]
                deriv.mu.hat.fin <- eval.fin[[6]]
                vcov.par.mu.hat <- vcov%*%t(deriv.mu.hat.fin)
                var.mu.hat <- diag(sigma.hat.fin^2)+(deriv.mu.hat.fin%*%vcov.par.mu.hat)
                mu.hat.df <- data.frame(cluster=clust,mu.hat=mu.hat.fin,sigma.hat=sigma.hat.fin)
                vcov.fix.mu.hat <- vcov.par.mu.hat[-which.rdm,,drop=FALSE]
            }
            else {
                mu.hat.df <- NA
                var.mu.hat <- NA
                vcov.fix.mu.hat <- NA
            }
        }
        else {
            param.fin <- rep(NA,n.par)
            code.fin <- 9
            loglik <- NA
            iterations <- 0
            names(param.fin) <- param.names
            gradient=rep(NA,n.par)
            vcov <- matrix(NA,n.par,n.par)
            colnames(vcov) <- rownames(vcov) <- param.names
            mu.hat.df <- NA
            var.mu.hat <- NA
            vcov.fix.mu.hat <- NA
        }
    }
    else if (mode=="eval"){
        eval.fin <- LL.Tot(init,mu.hat.LT=1)
        loglik <- eval.fin[[1]]
        gradient <- eval.fin[[2]]
        hessian <- eval.fin[[3]]
        param.fin <- init
        names(param.fin) <- param.names
        code.fin <- 9
        iterations <- 0
        vcov <- matrix(NA,n.par,n.par)
        colnames(vcov) <- rownames(vcov) <- param.names
        mu.hat.df <- NA
        var.mu.hat <- NA
        vcov.fix.mu.hat <- NA
        if (!is.null(random)) {
            mu.hat.df <- data.frame(cluster=clust,mu.hat=eval.fin[[4]],sigma.hat=eval.fin[[5]])
        }
    }
    time1 <- as.numeric(proc.time()[3])
    PrintData <- data.frame(Name=name.data,N.Obs.Tot=n.obs.tot,N.Obs=n.obs,N.Events=n.events,N.Clust=n.clust,row.names="")
    PrintDetails <- data.frame(Iter=iterations,Eval=parent.neval, Base=base,Nb.Leg=n.gleg,Nb.Aghq=n.aghq,Optim="nlm",Method="---",Code=code.fin,LogLik=loglik,Total.Time=(time1-time0),row.names="")

    ## Part of the results printed on screen
    cat("\nData\n")
    print(PrintData)
    cat("\nDetails\n")
    print(PrintDetails)

    res.FAR <- list(dataset=name.data,
         call=call,
         formula=formula,
         expected=ifelse(!is.null(expected),expected,NA),
         xlevels=Xlevels,
         n.obs.tot=n.obs.tot,
         n.obs=n.obs,
         n.events=n.events,
         n.clust=ifelse(!is.null(random),n.clust,1),
         n.time.0=ifelse(Survtype=="right",length(Idx.Time.0),0),
         base=base,
         max.time=max.time,
         boundary.knots=c(BoI,BoS),
         degree=degree[1],
         knots=knots,
         names.ph=names.fix,
         random=ifelse(!is.null(random),random,NA),
         init=init,
         coefficients=param.fin,
         std.errors=sqrt(diag(vcov)),
         vcov=vcov,
         gradient=gradient,
         hessian=hessian,
         mu.hat=mu.hat.df,
         var.mu.hat=var.mu.hat,
         vcov.fix.mu.hat=t(vcov.fix.mu.hat),
         n.par=n.par,
         n.gleg=n.gleg,
         n.aghq=n.aghq,
         fnoptim="nlm",
         method=NA,
         code=code.fin,
         loglik=loglik,
         iter=iterations,
         eval=parent.neval,
         time.elapsed=time1-time0)

    class(res.FAR) <- "mexhaz"
    res.FAR

}

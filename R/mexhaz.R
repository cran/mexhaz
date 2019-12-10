mexhaz <- function(formula,data,expected=NULL,base=c("weibull","exp.bs","exp.ns","pw.cst"),degree=3,knots=NULL,bound=NULL,n.gleg=20,init=NULL,random=NULL,n.aghq=10,fnoptim=c("nlm","optim"),verbose=0,method="Nelder-Mead",iterlim=10000,numHess=FALSE,print.level=1,...){

    time0 <- as.numeric(proc.time()[3])
    FALCenv <- environment()

    call <- match.call()
    base <- match.arg(base)
    fnoptim <- match.arg(fnoptim)

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

    if (base=="exp.bs" & !degree%in%c(1:3)){
        stop("This function can only be used to estimate log-hazards described by B-splines of degree 1 to 3...")
    }


    # Function that controls what is printed during the optimisation procedure
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
    if (n.gleg<=0 | round(n.gleg,0)!=n.gleg){
        stop("The 'n.gleg' argument must be a strictly positive integer.")
    }
    gl <- gauss.quad(n=n.gleg,kind="legendre")
    gln <- gl$nodes
    lglw <- log(gl$weights)

    # Data preparation
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
        if (base=="pw.cst"){
            HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HazardBs0R,x,nph,timecat,fixobs,param,paramf,as.double(matk))
            }
        }
        else if (base=="exp.bs" & degree==1){
            HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HazardBs1R,x,nph,timecat,fixobs,param,paramf,as.double(matk),as.double(totk))
            }
        }
        else if (base=="exp.bs" & degree%in%c(2,3)){
            HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HazardBs23R,x,nph,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
            }
        }
        else if (base=="exp.ns"){
            HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HazardNsR,x,nph,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2))
            }
        }
    }
    else if (ncol(Y)==3){
        if (Survtype!="counting"){
            stop(paste("mexhaz does not support \"", Survtype, "\" type of censoring with (time, time2] survival data",
                       sep = ""))
        }
        time.obs <- Y[,1:2]   # Entry time / Follow-up time
        status.obs <- Y[,3]   # Status variable
        if (is.null(random)){
            if (base=="pw.cst"){
                HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HazardBs0C,x0,x,nph,timecat0,timecat,fixobs,param,paramf,as.double(matk))
                }
            }
            else if (base=="exp.bs" & degree==1){
                HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HazardBs1C,x0,x,nph,timecat0,timecat,fixobs,param,paramf,as.double(matk),as.double(totk))
                }
            }
            else if (base=="exp.bs" & degree%in%c(2,3)){
                HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HazardBs23C,x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
                }
            }
            else if (base=="exp.ns"){
                HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HazardNsC,x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2))
                }
            }
        }
        else {
            if (base=="pw.cst"){
                HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HazardBs0L,x0,x,nph,timecat0,timecat,fixobs,param,paramf,as.double(matk))
                }
            }
            else if (base=="exp.bs" & degree==1){
                HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HazardBs1L,x0,x,nph,timecat0,timecat,fixobs,param,paramf,as.double(matk),as.double(totk))
                }
            }
            else if (base=="exp.bs" & degree%in%c(2,3)){
                HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HazardBs23L,x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
                }
            }
            else if (base=="exp.ns"){
                HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                    .Call(C_HazardNsL,x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2))
                }
            }
        }
    }

    if (base=="weibull"){
        if (!(Survtype=="counting" & !is.null(random))){
            HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                l.lambda.beta <- NULL
                Lambda.beta <- NULL
                lx <- log(x)
                log.p.LT.1 <- paramf[1] + as.vector(fixobs%*%paramf[-1])
                log.p.LT.2 <- param[1] + as.vector(nph%*%param[-1])
                l.lambda.beta <- log.p.LT.2+lx*(exp(log.p.LT.2)-1)+log.p.LT.1
                if (is.null(x0)){
                    Lambda.beta <- exp(log.p.LT.1)*x^exp(log.p.LT.2)
                }
                else {
                    Lambda.beta <- exp(log.p.LT.1)*(x^exp(log.p.LT.2)-x0^exp(log.p.LT.2))
                }
                valtot <- sum(l.lambda.beta) + sum(Lambda.beta)
                Test <- sum((is.nan(valtot)) | (valtot==Inf))
                Result <- list(LogHaz=l.lambda.beta, HazCum0=0, HazCum=Lambda.beta, Test=Test)
                return(Result)
            }
        }
        else {
            HazardInt <- function(x0,x,nph,timecat0,timecat,fixobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                l.lambda.beta <- NULL
                Lambda.beta <- NULL
                lx <- log(x)
                log.p.LT.1 <- paramf[1] + as.vector(fixobs%*%paramf[-1])
                log.p.LT.2 <- param[1] + as.vector(nph%*%param[-1])
                l.lambda.beta <- log.p.LT.2+lx*(exp(log.p.LT.2)-1)+log.p.LT.1
                Lambda.beta0 <- exp(log.p.LT.1)*x0^exp(log.p.LT.2)
                Lambda.beta  <- exp(log.p.LT.1)*x^exp(log.p.LT.2)
                valtot <- sum(l.lambda.beta) + sum(Lambda.beta)
                Test <- sum((is.nan(valtot)) | (valtot==Inf))
                Result <- list(LogHaz=l.lambda.beta, HazCum0=Lambda.beta0, HazCum=Lambda.beta, Test=Test)
                return(Result)
            }
        }
    }

    if (Survtype=="counting" & !is.null(random)){
        Frailty.Adapt <- function(nodes, nodessquare, logweights, clust, clustd, expect, betal, betaL0, betaL, A0, A, var, mh0, muhatcond){
            obj <- .Call(C_FrailtyAdaptL, nodes, nodessquare, logweights, clust, clustd, expect, betal, betaL0, betaL, A0, A, var, mh0, muhatcond)
        }
    }
    else {
        Frailty.Adapt <- function(nodes, nodessquare, logweights, clust, clustd, expect, betal, betaL0, betaL, A0, A, var, mh0, muhatcond){
            obj <- .Call(C_FrailtyAdapt, nodes, nodessquare, logweights, clust, clustd, expect, betal, betaL, A, var, muhatcond)
        }
    }

    data.fix <- model.frame(FormulaF,data=data,na.action=na.pass)
    data.nph <- model.frame(FormulaN,data=data,na.action=na.pass)
    n.obs.tot <- dim(data)[1] # Total number of observations

    # Expected hazard (a vector of 0 is created if NULL)
    if (!is.null(expected)){
        lambda.pop <- data[,expected]
        if (min(lambda.pop,na.rm=TRUE)<0){
            stop("The expected hazard for some observations is negative!")
        }
        withExp <- 1
    }
    else {
        lambda.pop <- rep(0,n.obs.tot)
        withExp <- 0
    }

    # Remove observations containing missing values and perform several formatting operations on the data
    if (!is.null(random)){
        random.obs <- data[,random]
        Idx.NA.Rdm <- which(is.na(random.obs))
        if (length(Idx.NA.Rdm)>0){
            warning("Cluster information was missing for some observations. These observations were consequently removed...")
            random.obs <- random.obs[-Idx.NA.Rdm]
            time.obs <- time.obs[-Idx.NA.Rdm,,drop=FALSE]
            status.obs <- status.obs[-Idx.NA.Rdm]
            lambda.pop <- lambda.pop[-Idx.NA.Rdm]
            data.fix <- data.fix[-Idx.NA.Rdm,,drop=FALSE]
            data.nph <- data.nph[-Idx.NA.Rdm,,drop=FALSE]
        }
        random.obs <- as.factor(random.obs)
        clust <- levels(random.obs)
        IdxRnd <- order(random.obs)
        # The dataset MUST be ordered by the levels of the clustering variable
        random.obs <- random.obs[IdxRnd]
        time.obs <- time.obs[IdxRnd,,drop=FALSE]
        status.obs <- status.obs[IdxRnd]
        lambda.pop <- lambda.pop[IdxRnd]
        data.fix <- data.fix[IdxRnd,,drop=FALSE]
        data.nph <- data.nph[IdxRnd,,drop=FALSE]
        listIdx <- sapply(levels(random.obs),function(x){which(random.obs==x)},simplify=FALSE)
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

    # Remove unnecessary objects
    rm(TotData,data.fix,data.nph)

    # Creation of the different objects necessary to compute the hazards

    # Weibull hazard
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
        fix.obs <- fix.obs[,-1,drop=FALSE]
        nph.obs <- nph.obs[,-1,drop=FALSE]
        param.names <- c("logLambda","logRho",names.fix,names.nph)
        n.par.fix <- n.td.base+n.ntd+n.td.nph
        param.init <- rep(0,n.par.fix)
        param.init[1:2] <- log(0.1)
    }

    # Hazard modelled by the exponential of a B-spline / restricted cubic B-spline / Piecewise constant
    else if (base%in%c("exp.bs","exp.ns","pw.cst")){
        # Baseline hazard-related objects
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

            # For non time-dependent effects
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

            # For non time-dependent effects
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

        # For time-dependent effects
        n.td.nph <- length(names.nph)*n.td.base
        if (n.td.nph>0){
            which.td <- c(n.inter+1:n.td.base,(n.td.base+n.ntd)+1:n.td.nph)
        }
        else {
            which.td <- c(n.inter+1:n.td.base)
        }
        names.nph <- unlist(sapply(names.nph,function(x){paste(x,names.base,sep="*")}))

        fix.obs <- t(fix.obs) # Matrix of fixed effects has to be transposed for use by the Int/Delta functions
        nph.obs <- t(nph.obs) # Matrix of time-dependent effects has to be transposed for use by the Int/Delta functions
        param.names <- c(intercept,names.base,names.fix,names.nph)
        n.par.fix <- n.td.base+n.ntd+n.td.nph
        param.init <- rep(0,n.par.fix)
        param.init[1:(n.td.base+n.inter)] <- -1
    }

    # Creation of objects related to the random effect
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
        lambda.pop.delta <- lambda.pop[status.one]
        n.clust <- length(clust) # Number of clusters
        n.by.clust <- as.vector(table(random.obs)) # Size of each cluster
        n.by.clust.delta <- as.vector(table(random.obs[status.one]))
        Idx.clust.no.obs <- which(n.by.clust==0)
        if (length(Idx.clust.no.obs)>0){
            warning("Some clusters had no non-missing values for some covariables. These observations were consequently removed...")
            cat("The following clusters had no non-missing values for some covariables:\n")
            print(names(table(random.obs))[Idx.clust.no.obs])
            n.by.clust <- n.by.clust[-Idx.clust.no.obs]
            n.by.clust.delta <- n.by.clust.delta[-Idx.clust.no.obs]
            clust <- clust[-Idx.clust.no.obs]
            n.clust <- length(clust)
        }
        parent.cst.adj <- rep(0,n.clust)
        if (Survtype=="counting"){
            parent.cst.adj0 <- rep(0,n.clust)
        }
        else {
            parent.cst.adj0 <- NULL
        }
        param.names <- c(param.names,paste(random," [log(sd)]",sep=""))

        # Creation of the points and weights of Gauss-Hermite quadrature
        if (n.aghq<=0 | round(n.aghq,0)!=n.aghq){
            stop("The 'n.aghq' argument must be a strictly positive integer.")
        }
        gq <- gauss.quad(n=n.aghq,kind="hermite")
        x.H <- gq$nodes
        x.H.2 <- x.H^2
        log.rho.H <- log(gq$weights)
        log.rho.H[log.rho.H==-Inf] <- -.Machine$double.xmax
    }


    # Initial values
    if (!is.null(init)){
        if (length(init)!=n.par){
            print(data.frame(Time.Dep.baseline=n.td.base,Non.TD.Effects=n.ntd,TD.Effects=n.td.nph,Random.Effect=n.rand))
            stop("Wrong number of init parameters",call.=FALSE)
        }
    }
    else {
        init <- c(param.init,rep(0.1,n.rand))
    }

    # Creation of variables that will be modified by the Hazard and LL.Tot functions
    parent.neval <- 0
    parent.param <- init
    parent.logLik <- -.Machine$double.xmax

    # Function that actually computes the log-likelihood
    LL.Tot <- function(p.LT,mu.hat.LT=0,env=FALCenv){

        p.td <- p.LT[which.td]
        p.ntd <- p.LT[which.ntd]

        temp.H <- HazardInt(x0=time.obs.0,x=time.obs,nph=nph.obs,timecat0=time.cat.0,timecat=time.cat,fixobs=fix.obs,param=p.td,paramf=p.ntd,deg=degree,n=gln,lw=lglw,matk=MatK,totk=vec.knots,intk=int.knots,nsadj1=NsAdj[[1]],nsadj2=NsAdj[[2]])
        if (temp.H$Test){
            res.LT <- .Machine$double.xmax
        }
        else {
            if (!is.null(random)){
                var.w <- exp(2*p.LT[which.rdm])
                if (Survtype=="counting"){
                    ## Mode of the integrand of the denominator of the cluster-specific marginal likelihood
                    MH0 <- -lambertW0(unlist(lapply(listIdx,function(x){sum(temp.H$HazCum0[x])*var.w})))
                }
                else {
                    MH0 <- 0
                }
                temp.LT <- Frailty.Adapt(nodes=x.H, nodessquare=x.H.2, logweights=log.rho.H, clust=n.by.clust, clustd=n.by.clust.delta, expect=lambda.pop.delta, betal=temp.H$LogHaz[status.one], betaL0=temp.H$HazCum0, betaL=temp.H$HazCum, A0=parent.cst.adj0, A=parent.cst.adj, var=var.w, mh0=MH0, muhatcond=mu.hat.LT)
                if (mu.hat.LT==1)
                    res.LT <- temp.LT$MuHat
                else if (mu.hat.LT==2)
                    res.LT <- temp.LT$SigmaHat
                else {
                    env$parent.cst.adj0 <- temp.LT$CstAdj0
                    env$parent.cst.adj <- temp.LT$CstAdj
                    res.LT <- temp.LT$LogLik
                    res.LT[is.nan(res.LT) | abs(res.LT)==Inf] <- .Machine$double.xmax
                }
            }
            else {
                if (withExp==1){
                    temp.l <- exp(temp.H$LogHaz)+lambda.pop
                    if (sum(temp.l-lambda.pop)==0){
                        res.LT <- .Machine$double.xmax
                    } # Prevents a kind of catastrophic cancellation
                    else {
                        log.lambda <- log(temp.l)
                        log.lambda[log.lambda==Inf] <- .Machine$double.xmax
                        res.LT <- sum(temp.H$HazCum - status.obs*log.lambda)
                        res.LT <- min(res.LT,.Machine$double.xmax)
                    }
                }
                else {
                    log.lambda <- temp.H$LogHaz
                    log.lambda[log.lambda==Inf] <- .Machine$double.xmax
                    res.LT <- sum(temp.H$HazCum - status.obs*log.lambda)
                    res.LT <- min(res.LT,.Machine$double.xmax)
                }
            }
            if (mu.hat.LT==0){
                if (sum(p.LT-parent.param)){
                    verbose.ll()
                    env$parent.neval <- parent.neval + 1
                }
                env$parent.param <- p.LT + 0 # Necessary...
                env$parent.logLik <- res.LT
            }
        }
        return(res.LT)
    }

    # Launching the optimisation procedure
    if (fnoptim=="nlm"){
        time0 <- as.numeric(proc.time()[3])
        mod.lik <- nlm(LL.Tot,init,iterlim=iterlim,hessian=TRUE,print.level=print.level,...)
        param.fin <- mod.lik$estimate
        code.fin <- mod.lik$code
        loglik <- -mod.lik$minimum
        iterations <- mod.lik$iterations
    }
    else if (fnoptim=="optim"){
        time0 <- as.numeric(proc.time()[3])
        mod.lik <- optim(init,LL.Tot,method=method,hessian=TRUE,...)
        param.fin <- mod.lik$par
        code.fin <- mod.lik$convergence
        loglik <- -mod.lik$value
        iterations <- "---"
    }
    names(param.fin) <- param.names

    # Hessian matrix
    verbose.ll <- function(){} # Suppresses iteration printing
    cat("Computation of the Hessian\n")
    if (numHess){
        hessian.fin <- hessian(LL.Tot,param.fin,mu.hat.LT=0,method.args=list(r=4))
    }
    else {
        hessian.fin <- mod.lik$hessian
    }
    vcov <- try(solve(hessian.fin),silent=TRUE)
    if (is.character(vcov)){
        vcov <- matrix(NA,n.par,n.par)
        warning("Unable to invert the Hessian matrix...")
    }
    colnames(vcov) <- rownames(vcov) <- param.names

    # Empirical Bayes Estimates
    if (!is.null(random)) {
        cat("Computation of the covariance matrix of the shrinkage estimators\n")
        mu.hat.fin <- LL.Tot(param.fin,mu.hat.LT=1)
        sigma.hat.fin <- LL.Tot(param.fin,mu.hat.LT=2)
        deriv.mu.hat.fin <- jacobian(LL.Tot,param.fin,mu.hat.LT=1,method.args=list(r=4))
        vcov.par.mu.hat <- vcov%*%t(deriv.mu.hat.fin)
        var.mu.hat <- diag(sigma.hat.fin^2)+(deriv.mu.hat.fin%*%vcov.par.mu.hat)
        mu.hat.df <- data.frame(cluster=clust,mu.hat=mu.hat.fin)
        vcov.fix.mu.hat <- vcov.par.mu.hat[-which.rdm,,drop=FALSE]
    }
    else {
        mu.hat.df <- 0
        var.mu.hat <- 0
        vcov.fix.mu.hat <- 0
    }

    time1 <- as.numeric(proc.time()[3])
    PrintData <- data.frame(Name=name.data,N.Obs.Tot=n.obs.tot,N.Obs=n.obs,N.Events=n.events,N.Clust=n.clust,row.names="")
    PrintDetails <- data.frame(Iter=iterations,Eval=parent.neval+1,Base=base,Nb.Leg=n.gleg,
               Nb.Aghq=n.aghq,Optim=fnoptim,Method=ifelse(fnoptim=="optim",method,"---"),
               Code=code.fin,LogLik=loglik,Total.Time=(time1-time0),row.names="")

    # Part of the results printed on screen
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
         coefficients=param.fin,
         std.errors=sqrt(diag(vcov)),
         vcov=vcov,
         mu.hat=mu.hat.df,
         var.mu.hat=var.mu.hat,
         vcov.fix.mu.hat=t(vcov.fix.mu.hat),
         n.par=n.par,
         n.gleg=n.gleg,
         n.aghq=n.aghq,
         fnoptim=fnoptim,
         method=ifelse(fnoptim=="optim",method,NA),
         code=code.fin,
         loglik=loglik,
         iter=iterations,
         eval=parent.neval+1,
         time.elapsed=time1-time0)

    class(res.FAR) <- "mexhaz"
    res.FAR

}

residuals.mexhaz <- function(object,data,hazExpect=NULL,cumExpect=NULL,...){

    time0 <- as.numeric(proc.time()[3])
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
    excess <- object$expected

    if (!is.na(object$random)){
        stop("Currently, the 'residuals.mexhaz' function is not implemented for a model including a random effect...")
    }

    data <- data[!1:object$n.obs.tot%in%object$idx.NA,]
    n.obs <- dim(data)[1]
    if (n.obs!=object$n.obs){
        stop("Something went wrong when trying to recreate the dataset. Please, try to clean the dataset before fitting the model...")
    }

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

    FormulaM <- update(tot.formula,".~NULL")
    FormulaF <- update(tot.formula,paste("NULL~.",nTerm2,sep="-"))
    Xlevel.formula <- terms(FormulaF,data=data)
    m <- eval(model.frame(FormulaM,data=data))
    Y <- model.extract(m,"response")
    if (!inherits(Y,"Surv")){
        stop("Response must be a Surv() object...")
    }
    Survtype <- attr(Y, "type")
    if (ncol(Y)==3){
        stop("residuals.mexhaz does not support data with delayed entry for the time being...")
    }
    if (ncol(Y)==2){
        if (Survtype!="right"){
            stop(paste("residuals.mexhaz does not support \"", Survtype, "\" type of censoring with (0, time] survival data",
                       sep = ""))
        }
        time.pts <- Y[,1,drop=FALSE]   # Follow-up time
        event <- Y[,2]   # Status variable
    }

    time.pts[time.pts==0] <- 1/730.5

    if (!is.na(excess)){
        if (!(!is.null(hazExpect) & !is.null(cumExpect))){
            stop("Calculation of residuals for excess hazard models requires that both the expected hazard and the expected cumulative hazard be provided in the dataset (their names being specified through the arguments 'hazExpect' and 'cumExpect')...")
        }
        else {
            expect <- data[,hazExpect]
            cumexp <- data[,cumExpect]
        }
    }
    else {
        expect <- cumexp <- rep(0,dim(data)[1])
    }

    FALCenv <- environment()

    if (is.na(excess)){
        if (base=="pw.cst"){
            HazardInt <- function(x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HGH_PwR,x,nph,timecat,fixobs,param,paramf,as.double(matk))
            }
        }
        else if (base=="exp.bs" & degree==1){
            HazardInt <- function(x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HGH_BsR,x,nph,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
            }
        }
        else if (base=="exp.bs" & degree%in%c(2,3)){
            HazardInt <- function(x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HGH_BsR,x,nph,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
            }
        }
        else if (base=="exp.ns"){
            HazardInt <- function(x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HGH_NsR,x,nph,timecat,fixobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2))
            }
        }
        else if (base=="weibull"){
            HazardInt <- function(x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HGH_WbR,x,nph,fixobs,param,paramf)
            }
        }
    }
    else {
        if (base=="pw.cst"){
            HazardInt <- function(x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HGH_PwRx,x,nph,timecat,fixobs,lambdaobs,param,paramf,as.double(matk))
            }
        }
        else if (base=="exp.bs" & degree==1){
            HazardInt <- function(x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HGH_BsRx,x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
            }
        }
        else if (base=="exp.bs" & degree%in%c(2,3)){
            HazardInt <- function(x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HGH_BsRx,x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk))
            }
        }
        else if (base=="exp.ns"){
            HazardInt <- function(x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HGH_NsRx,x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,as.double(matk),as.double(totk),as.double(intk),as.double(nsadj1),as.double(nsadj2))
            }
        }
        else if (base=="weibull"){
            HazardInt <- function(x,nph,timecat,fixobs,lambdaobs,param,paramf,deg,n,lw,matk,totk,intk,nsadj1,nsadj2){
                .Call(C_HGH_WbRx,x,nph,fixobs,lambdaobs,param,paramf)
            }
        }
    }

    Reorder <- function(vec,mat,idx){
        mat[idx] <- vec
        mat <- t(mat)
        mat[idx] <- vec
        return(mat)
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
    gl <- gauss.quad(n=n.gleg,kind="legendre")
    gln <- gl$nodes
    lglw <- log(gl$weights)

    ## Formatting the data
    fix.new <- model.matrix(FormulaF,data=data)
    names.fix <- colnames(fix.new)[-1]
    nph.new <- model.matrix(FormulaN,data=data)
    nbtd <- dim(nph.new)[2]
    names.nph <- colnames(nph.new)[-1]

    ## Creation of the different objects necessary to compute the hazards

    ## Weibull hazard
    if (base=="weibull"){
        time.cat <- NULL
        vec.knots <- NULL
        int.knots <- NULL
        MatK <- NULL
        NsAdj <- list(NULL,NULL)
        time.new <- time.pts
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
        fix.new.2 <- fix.new[,-1,drop=FALSE]
        nph.new.2 <- nph.new[,-1,drop=FALSE]
    }

    ## Hazard modelled by the exponential of a B-spline
    else if (base%in%c("exp.bs","exp.ns","pw.cst")){

        cuts <- c(0,knots,max.time)
        time.cat.lab <- cut(time.pts,breaks=cuts,include.lowest=TRUE)
        time.cat <- as.numeric(time.cat.lab)-1

        if (base=="pw.cst"){
            time.new <- time.pts-cuts[time.cat+1]
            vec.knots <- NULL
            int.knots <- NULL
            MatK <- c(knots,BoS)-c(0,knots)
            NsAdj <- list(NULL,NULL)
            n.td.base <- length(knots)+1

            ## For non time-dependent effects
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
                time.cat.lab <- cut(time.pts,breaks=int.knots,include.lowest=TRUE)
                time.cat <- as.numeric(time.cat.lab)-1
                degree <- c(degree,ltk,firstK) # ltk and firstK are passed to the HazardInt function through the 'degree' argument which is not used with NS
            }
            n.td.base <- dbase + length(knots)

            ## For non time-dependent effects
            n.ntd <- dim(fix.new)[2]
            if (n.ntd>1){
                which.ntd <- c(1,(n.td.base+1)+1:(n.ntd-1))
            }
            else {
                which.ntd <- 1
            }
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

    }

    fix.new <- t(fix.new) # Matrix of prop effects has to be transposed for use by the Int/Delta functions
    nph.new <- t(nph.new) # Matrix of time-dep effects has to be transposed for use by the Int/Delta functions


    OrderObs <- order(time.pts)
    time.ord <- time.pts[OrderObs]
    time.nord <- time.new[OrderObs]
    time.cord <- time.cat[OrderObs]
    unique.time <- sort(unique(time.ord))
    n.ut <- length(unique.time)
    cut.T <- cut(time.pts,breaks=c(unique.time[1]-1,unique.time,unique.time[n.ut]+1))
    idx.obs <- as.integer(cut.T)
    idx.time <- which(diff(c(0,idx.obs[OrderObs]))==1)
    nb.par <- length(c(which.ntd,which.td))
    temp.H <- HazardInt(x=time.new,nph=nph.new,timecat=time.cat,fixobs=fix.new,lambdaobs=expect,param=coef[which.td],paramf=coef[which.ntd],deg=degree,n=gln,lw=lglw,matk=MatK,totk=vec.knots,intk=int.knots,nsadj1=NsAdj[[1]],nsadj2=NsAdj[[2]])
    temp.H0 <- HazardInt(x=time.nord[idx.time],nph=rep(1,n.ut),timecat=time.cord[idx.time],fixobs=rep(1,n.ut),lambdaobs=rep(0,n.ut),param=coef[which.td],paramf=coef[which.ntd[1]],deg=degree,n=gln,lw=lglw,matk=MatK,totk=vec.knots,intk=int.knots,nsadj1=NsAdj[[1]],nsadj2=NsAdj[[2]])
    lambda <- exp(temp.H$LogHaz)
    lambda0 <- exp(temp.H0$LogHaz)
    Order <- order(c(which.ntd,which.td))

    dlhi.db <- matrix(temp.H$GradLogHaz,n.obs,nb.par)
    dHi.db <- matrix(temp.H$GradCum,n.obs,nb.par)
    d2lhi.db2 <- event*matrix(temp.H$HessLHaz,n.obs,0.5*nb.par*(nb.par+1))
    d2Hi.db2 <- matrix(temp.H$HessCum,n.obs,0.5*nb.par*(nb.par+1))
    ScoreIa <- dHi.db[,Order]
    ScoreIb <- (event*dlhi.db)[,Order]
    ScoreI <- (event*dlhi.db-dHi.db)[,Order]
    Score <- (event%*%dlhi.db-apply(dHi.db,2,sum))[Order]
    Mat <- matrix(0,nb.par,nb.par)
    IdxM <- which(lower.tri(Mat,diag=TRUE))
    Hess <- Reorder(apply(d2Hi.db2-d2lhi.db2,2,sum),Mat,IdxM)[Order,Order]

    res0 <- data.frame(time=unique.time,hazard0=lambda0,cumhaz0=temp.H0$HazCum)
    gradhaz0 <- matrix(temp.H0$GradLogHaz*lambda0,n.ut,1+length(which.td))
    gradcum0 <- matrix(temp.H0$GradCum,n.ut,1+length(which.td))
    res1 <- data.frame(time=time.pts,idxobs=idx.obs,event=event,hazard=lambda,cumhaz=temp.H$HazCum,hazexp=expect,cumexp=cumexp)
    covariates <- (t(fix.new))[,-1,drop=FALSE]
    if (length(covariates)==0){
        covariates <- NA
        exp.bx <- 1
    }
    else {
        exp.bx <- exp(covariates%*%coef[which.ntd[-1]])
    }

    res.RS <- list(call=call,
                   order.time=OrderObs,
                   results0=res0,
                   gradhaz0=gradhaz0,
                   gradcum0=gradcum0,
                   results=res1,
                   exp.bx=exp.bx,
                   covariates=covariates,
                   martingale.resid=event-temp.H$HazCum-cumexp,
                   score.resid=ScoreI,
                   scoreR1=ScoreIa,
                   scoreR2=ScoreIb,
                   score=Score,
                   hessian=Hess)
    res.RS <- res.RS[!is.na(res.RS)]

    class(res.RS) <- "residMexhaz"
    res.RS
}

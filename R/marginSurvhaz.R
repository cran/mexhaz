marginSurvhaz <- function(lA,B,C,dlA,dC,s2,n.aghq=100,grad=FALSE){

    gq <- gauss.quad(n=n.aghq,kind="hermite")
    x.H <- gq$nodes
    rho.H <- gq$weights
    Bs2 <- B*s2
    B2s2 <- B^2*s2
    Nu <- lambertW0(C*s2*exp(Bs2))
    Mu <- Bs2-Nu
    SumB <- function(x){
        temp <- sqrt(2*s2/(1+x))*x.H
        rho.H%*%exp(-x/s2*(exp(temp)-0.5*(1+(temp+1)^2)))
    }
    SumHB <- function(x){
        temp <- sqrt(2*s2/(1+x))*x.H
        g <- exp(temp)-0.5*(1+(temp+1)^2)
        gp <- x.H*(exp(temp)-(temp+1))
        gp2 <- x.H^2*(exp(temp)-1)
        exp.tg <- exp(-x/s2*g)
        rho.H%*%cbind(g*exp.tg,gp*exp.tg,gp2*exp.tg,g^2*exp.tg,g*gp*exp.tg,gp^2*exp.tg)
    }
    SnB <- sapply(Nu,SumB)
    Int <- exp(lA-0.5*log(pi*(1+Nu))-0.5*(Nu*(Nu+2)/s2-B2s2)+log(SnB))

    ## Gradient
    if (grad==TRUE){
        dNu.db <- dC*Nu/((1+Nu)*C)
        dNu.dls <- 2*Nu*(1+Bs2)/(1+Nu)

        Psi <- sqrt(2*s2/(1+Nu))
        dPsi.db <- -Psi*dNu.db/(2*(1+Nu))
        dPsi.dls <- -Psi*(dNu.dls/(2*(1+Nu))-1)
        dPsi.di <- cbind(dPsi.db,dPsi.dls)

        Theta <- -Nu/s2
        dTheta.di <- -(1/s2)*cbind(dNu.db,dNu.dls-2*Nu)

        C1 <- (1/(1+Nu)+2*(Nu+1)/s2)
        dF.db <- dlA-0.5*C1*dNu.db
        dF.dls <- -0.5*C1*dNu.dls+(Nu*(Nu+2)/s2+B2s2)
        dF.di <- cbind(dF.db,dF.dls)
        sumHB <- t(sapply(Nu,SumHB))/SnB
        dlSn.di <- dTheta.di*sumHB[,1]+Theta*dPsi.di*sumHB[,2]
        dInt.di <- Int*(dF.di+dlSn.di)
    }
    else {
        dInt.di <- NA
    }

    return(list(Int,dInt.di))
}

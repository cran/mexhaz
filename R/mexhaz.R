mexhaz <- function(formula, data, expected=NULL, base=c("weibull","exp.bs","exp.ns","pw.cst"), degree=3, knots=NULL, bound=NULL, n.gleg=20, init=NULL, random=NULL, n.aghq=10, recurrent=FALSE, fnoptim=c("nlm","optim"), verbose=0, method="Nelder-Mead", iterlim=10000, numHess=FALSE, print.level=1, exactGradHess=TRUE, gradtol=ifelse(exactGradHess,1e-8,1e-6), testInit=TRUE, keep.data=FALSE, ...){

    base <- match.arg(base)
    fnoptim <- match.arg(fnoptim)
    name.data <- paste0(substitute(data))
    if (length(name.data)>1){
        name.data <- name.data[2]
    }

    if (exactGradHess==TRUE & !(!is.null(expected) & !is.null(random))){
        res <- mexhazEgh(formula=formula,data=data,expected=expected,base=base,degree=degree,knots=knots,bound=bound,n.gleg=n.gleg,init=init,random=random,n.aghq=n.aghq,recurrent=recurrent,verbose=verbose,iterlim=iterlim,print.level=print.level,gradtol=gradtol,testInit=testInit,keep.data=keep.data,name.data=name.data,...)
    }
    else {
        res <- mexhazStd(formula=formula,data=data,expected=expected,base=base,degree=degree,knots=knots,bound=bound,n.gleg=n.gleg,init=init,random=random,n.aghq=n.aghq,fnoptim=fnoptim,verbose=verbose,method=method,iterlim=iterlim,numHess=numHess,print.level=print.level,gradtol=gradtol,keep.data=keep.data,name.data=name.data,...)
    }

    res$call <- match.call()
    return(res)
}

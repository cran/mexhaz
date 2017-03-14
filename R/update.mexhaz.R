update.mexhaz <- function (object, formula, data, expected = NULL, base = c("weibull",
    "exp.bs", "exp.ns", "pw.cst"), degree = 3, knots = NULL,
    bound = NULL, n.gleg = 20, init = NULL, random = NULL, n.aghq = 10,
    fnoptim = c("nlm", "optim"), verbose = 100, method = "Nelder-Mead",
    iterlim = 10000, numHess = FALSE, print.level = 1,...) {
    Call <- match.call()
    call.mod <- object$call
    if (!("formula"%in%names(Call))){
        Call$formula <- .~.
    }
    call.mod$formula <- update(as.formula(call.mod$formula),as.formula(Call$formula))
    new.call <- call.mod[is.na(match(names(call.mod),c("init")))]
    if (length(which(is.na(match(names(Call),c("","object","formula")))))>0){
        call.up <- Call[is.na(match(names(Call),c("","object","formula")))]
        if (length(which(!is.na(match(names(new.call),names(call.up)))))>0){
            new.call[!is.na(match(names(new.call),names(call.up)))] <- call.up[!is.na(match(names(call.up),names(new.call)))]
        }
        if (length(which(is.na(match(names(call.up),names(new.call)))))>0){
            new.call[names(call.up[is.na(match(names(call.up),names(new.call)))])] <- call.up[is.na(match(names(call.up),names(new.call)))]
        }
    }
    if (!("init"%in%names(new.call))){
        base <- eval(new.call$base)
        data <- eval(new.call$data)
        formula <- eval(new.call$formula)
        if (!("degree"%in%names(new.call))){
            degree <- object$degree
        }
        else degree <- eval(new.call$degree)
        if (!("knots"%in%names(new.call))){
            knots <- object$knots
        }
        else knots <- eval(new.call$knots)
        if (!("random"%in%names(new.call))){
            random <- NULL
        }
        else random <- eval(new.call$random)
        if (base=="exp.bs"){
            dbase <- degree
        }
        if (base=="exp.ns"){
            degree <- 3
            dbase <- 1
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
        FormulaF <- update(tot.formula,paste("NULL~.",nTerm2,sep="-"))
        data.fix <- model.frame(FormulaF,data=data,na.action=na.pass)
        data.nph <- model.frame(FormulaN,data=data,na.action=na.pass)
        fix.obs <- model.matrix(FormulaF,data=data.fix,drop.unused.levels=TRUE)
        names.fix <- colnames(fix.obs)[-1]
        nph.obs <- model.matrix(FormulaN,data=data.nph,drop.unused.levels=TRUE)
        nbtd <- dim(nph.obs)[2]
        names.nph <- colnames(nph.obs)[-1]
        # Weibull hazard
        if (base=="weibull"){
            degree <- NA
            n.td.base <- 1
            n.ntd <- dim(fix.obs)[2]
            n.td.nph <- length(names.nph)
            if (n.td.nph>0){
                names.nph <- paste("Rho",names.nph,sep="*")
            }
            param.names <- c("Lambda","Rho",names.fix,names.nph)
            n.par.fix <- n.td.base+n.ntd+n.td.nph
            param.init <- rep(0,n.par.fix)
            param.init[1:2] <- 0.1
        }

        # Hazard modelled by the exponential of a B-spline / restricted cubic B-spline / Piecewise constant
        else if (base%in%c("exp.bs","exp.ns","pw.cst")){
            # Baseline hazard-related objects
            if (!is.null(knots)){
                if (sum(abs(knots-knots[order(knots)]))>0){
                    knots <- knots[order(knots)]
                }
                if (length(unique(knots))<length(knots)){
                    knots <- unique(knots)
                }
            }
            if (base=="pw.cst"){
                degree <- 0
                cuts <- c(0,knots,object$max.time)
                n.td.base <- length(knots)+1
                names.base <- levels(cut(0,breaks=cuts,include.lowest=TRUE))
                # For non time-dependent effects
                fix.obs <- fix.obs[,-which(colnames(fix.obs)%in%colnames(nph.obs)),drop=FALSE]
                names.fix <- colnames(fix.obs)
                if (!is.null(names.fix)){
                    n.ntd <- dim(fix.obs)[2]
                }
                else {
                    n.ntd <- 0
                }
                intercept <- NULL
                n.inter <- 0
            }
            else if (base%in%c("exp.bs","exp.ns")){
                n.td.base <- dbase + length(knots)
                names.base <- paste(ifelse(base=="exp.bs","BS","NS"),degree,".",1:n.td.base,sep = "")
                # For non time-dependent effects
                n.ntd <- dim(fix.obs)[2]
                intercept <- "Intercept"
                n.inter <- 1
            }
            # For time-dependent effects
            n.td.nph <- length(names.nph)*n.td.base
            names.nph <- unlist(sapply(names.nph,function(x){paste(x,names.base,sep="*")}))
            param.names <- c(intercept,names.base,names.fix,names.nph)
            n.par.fix <- n.td.base+n.ntd+n.td.nph
            param.init <- rep(0,n.par.fix)
            param.init[1:(n.td.base+n.inter)] <- -1
        }

        n.rand <- 0
        if (!is.null(random)){
            n.rand <- 1
            param.names <- c(param.names,paste(random," (sd)",sep=""))
        }
        init <- c(param.init,rep(0.1,n.rand))
        if (length(which(!is.na(match(param.names,names(object$coefficients)))))>0){
            init[!is.na(match(param.names,names(object$coefficients)))] <- object$coefficients[!is.na(match(names(object$coefficients),param.names))]
        }
        new.call$init <- init
    }
    eval(new.call)
}

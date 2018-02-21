### sCorrect.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2018 (14:29) 
## Version: 
## Last-Updated: feb 21 2018 (18:15) 
##           By: Brice Ozenne
##     Update #: 628
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - sCorrect
#' @title  Compute the Derivative of the Information Matrix
#' @description Compute the derivative of the information matrix.
#' @name sCorrect
#'
#' @param object,x a \code{gls}, \code{lme}, or \code{lvm} object.
#' @param cluster [vector] the grouping variable relative to which the observations are iid.
#' Only required for \code{gls} models with no correlation argument.
#' @param vcov.param [matrix] the variance-covariance matrix of the estimates.
#' @param bias.correct,value [logical] should the standard errors of the coefficients be corrected for small sample bias?
#' @param numeric.derivative [logical] should a numerical derivative be used to compute the first derivative of the information matrix?
#' Otherwise an analytic formula is used.
#' @param return.score [internal] export the score.
#' @param ... [internal] only used by the generic method or by the <- methods.
#'
#' @concept small sample inference
#' @concept derivative of the score equation
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#' 
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- sequence.lvm(0)
#' set.seed(10)
#' d <- sim(m,n)
#'
#' ## linear model
#' e.lm <- lm(formula.lvm,data=d)
#' system.time(
#' sCorrect(e.lm) <- TRUE
#')
#' 
#' ## gls model
#' library(nlme)
#' e.gls <- gls(formula.lvm, data = d, method = "ML")
#' sCorrect(e.gls, cluster = 1:NROW(d)) <- TRUE
#' 
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' iid.tempo <- iid2(e.lvm, adjust.residuals = FALSE)
#' range(iid.tempo-iid(e.lvm))
#' ## difference due to the use of the observed info matrix vs. the expected one.
#' @export
`sCorrect` <-
  function(object, ...) UseMethod("sCorrect")


## * sCorrect.lm
#' @rdname sCorrect
#' @export
sCorrect.lm <- function(object, adjust.Omega = TRUE, adjust.n = TRUE,
                        score = TRUE, df = TRUE, numeric.derivative = FALSE,
                        param = NULL, data = NULL,
                        tol = 1e-5, n.iter = 20, ...){
    
    ## ** Extract quantities from object
    name.endogenous <- all.vars(stats::update(formula(object), ".~1"))
    n.cluster <- nobs(object) + length(object$na.action)
    
    if(is.null(param)){
        param <- .coef2(object)
        model.param <- param
    }else{
        model.param <- .coef2(object)
        if(any(names(param) %in% names(model.param) == FALSE)){
            stop("Argument \'param\' have appropriate names: \"",
                 paste(setdiff(names(param),names(model.param)), collapse = "\" \""),
                 "\" \n")
        }
        model.param[names(param)] <- param
    }
    if(is.null(data)){
        if(is.null(object$na.action)){
            object.residuals <- cbind(residuals(object))
            index.Omega <- NULL
        }else{
            object.residuals <- matrix(NA, nrow = n.cluster, ncol = 1,
                                       dimnames = list(NULL, name.endogenous))
            index.Omega <- setdiff(1:n.cluster, object$na.action)
            object.residuals[index.Omega,]  <- residuals(object)            
        }
    }else{
        object.residuals <- cbind(predict(object, newdata = data) - data[[name.endogenous]])
    }
    colnames(object.residuals) <- name.endogenous

    
    name.param <- names(model.param)
    name.meanparam <- attr(model.param,"mean.coef")
    name.varparam <- attr(model.param,"var.coef")
    object.sigma2 <- matrix(mean(object.residuals^2), nrow = 1, ncol = 1,
                            dimnames = list(name.endogenous, name.endogenous))
    
    ## ** Compute conditional moments
    dMoments <- conditionalMoment(object, name.endogenous = name.endogenous,
                                  second.order = df)

    ## ** correction
    if(df == FALSE){
        derivative <- "none"
    }else if(numeric.derivative){
        derivative <- "numeric"
    }else{
        derivative <- "analytic"
    }

    out <- .sCorrect(object,
                     param = model.param,
                     epsilon = object.residuals,
                     Omega = object.sigma2,
                     dmu = dMoments$dmu,
                     dOmega = dMoments$dOmega,
                     d2mu = dMoments$d2mu,
                     d2Omega = dMoments$d2Omega,
                     name.param = name.param,
                     name.meanparam = name.meanparam,
                     name.varparam = name.varparam,
                     name.endogenous = name.endogenous,
                     name.3deriv = dMoments$name.3deriv,
                     n.cluster = n.cluster,
                     index.Omega = index.Omega,
                     adjust.Omega = adjust.Omega,
                     adjust.n = adjust.n,
                     tol = tol,
                     n.iter = n.iter,
                     score = score,
                     derivative = derivative,
                     ...)
    
    ## ** export
    out$args <- list(adjust.Omega = adjust.Omega,
                     adjust.n = adjust.n,
                     df = df,
                     numeric.derivative = numeric.derivative,
                     tol = tol, n.iter = n.iter)
    return(out)    
}

## * sCorrect.gls
#' @rdname sCorrect
#' @export
sCorrect.gls <- function(object, cluster, adjust.Omega = TRUE, adjust.n = TRUE,
                         score = TRUE, df = TRUE, numeric.derivative = FALSE, 
                         param = NULL, data = NULL,
                         tol = 1e-5, n.iter = 20,
                         ...){

    ### ** limitations
    if(object$method!="ML"){
        warning("Implemented for Maximum Likelihood estimation not for REML\n")
    }
    ## check valid class for corStruct and varStruct: see .getVarCov2
    
    ### ** Extract quantities from the model
    ## *** data
    if(is.null(data)){
        data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE)
    }
    
    ## *** endogenous variable
    formula.object <- .getFormula2(object)
    name.Y <- all.vars(stats::update(formula.object, ".~1"))
    Y <- data[[name.Y]]
    
    ## *** parameters
    model.param <- .coef2(object)
    if(!is.null(param)){        
        if(any(names(param) %in% names(model.param) == FALSE)){
            stop("Argument \'param\' have appropriate names: \"",
                 paste(setdiff(names(param),names(model.param)), collapse = "\" \""),
                 "\" \n")
        }
        model.param[names(param)] <- param
    }
    name.param <- names(model.param)
    name.meanparam <- attr(model.param,"mean.coef")
    name.varparam <- c(attr(model.param,"var.coef"), attr(model.param,"cor.coef"))

    ## *** group
    res.cluster <- .getCluster2(object,
                                data = data,
                                cluster = cluster)

    ## *** repetition relative to each observation
    res.index <- .getIndexOmega2(object,
                                 param = model.param,
                                 attr.param = attributes(model.param),
                                 name.Y = name.Y,
                                 cluster = res.cluster$cluster,
                                 data = data)
    index.Omega <- res.index$index.Omega
    
    ### ** Reconstruct variance covariance matrix (residuals)
    Omega <- .getVarCov2.gls(object,
                             param = model.param,
                             attr.param = attributes(model.param),
                             name.endogenous = res.index$name.endogenous,
                             n.endogenous = res.index$n.endogenous,
                             ref.group = res.index$ref.group)
    
    cluster <- res.cluster$cluster
    n.cluster <- res.cluster$n.cluster
    name.endogenous <- res.index$name.endogenous
    n.endogenous <- res.index$n.endogenous

    ### ** Compute conditional moments and derivatives
    dMoments <- conditionalMoment(object,
                                  data = data,
                                  formula = formula.object,
                                  param = model.param,
                                  attr.param = attributes(model.param)[-1],
                                  ref.group = res.index$ref.group,
                                  second.order = df,
                                  n.cluster = n.cluster,
                                  cluster = cluster,
                                  name.endogenous = name.endogenous,
                                  index.Omega = index.Omega)

    ### ** Compute observed residuals
    epsilon <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                      dimnames = list(NULL, name.endogenous))
    for(iC in 1:n.cluster){
        iIndex <- which(cluster == iC)
        epsilon[iC,index.Omega[[iC]]] <- as.double(Y[iIndex] - dMoments$X[iIndex,,drop=FALSE] %*% model.param[name.meanparam])
    }    
    ## stats::residuals(object)-as.double(t(epsilon))

    ## ** Check missing value
    if(all(!is.na(epsilon))){
        index.Omega <- NULL
    }
    
    ## ** param with non-zero third derivative
    name.3deriv <- name.varparam
   
    ## ** correction
    if(df == FALSE){
        derivative <- "none"
    }else if(numeric.derivative){
        derivative <- "numeric"
    }else{
        derivative <- "analytic"
    }

    out <- .sCorrect(object,
                     param = model.param,
                     epsilon = epsilon,
                     Omega = Omega,
                     dmu = dMoments$dmu,
                     dOmega = dMoments$dOmega,
                     d2mu = NULL,
                     d2Omega = dMoments$d2Omega,
                     name.param = name.param,
                     name.meanparam = name.meanparam,
                     name.varparam = name.varparam,
                     name.endogenous = name.endogenous,
                     name.3deriv = name.3deriv,
                     n.cluster = n.cluster,
                     index.Omega = index.Omega,
                     adjust.Omega = adjust.Omega,
                     adjust.n = adjust.n,
                     tol = tol,
                     n.iter = n.iter,
                     score = score,
                     derivative = derivative,
                     ...)
    
    ## ** export
    out$args <- list(adjust.Omega = adjust.Omega,
                     adjust.n = adjust.n,                     
                     df = df,
                     cluster = cluster,
                     numeric.derivative = numeric.derivative,
                     tol = tol, n.iter = n.iter)
    return(out)          
 
}

## * sCorrect.lme
#' @rdname sCorrect
#' @export
sCorrect.lme <- sCorrect.gls

## * sCorrect.lvmfit
#' @rdname sCorrect
#' @export
sCorrect.lvmfit <- function(object, adjust.Omega = TRUE, adjust.n = TRUE,
                            score = TRUE, df = TRUE, numeric.derivative = FALSE, 
                            param = NULL, data = NULL,
                            tol = 1e-5, n.iter = 20,
                            ...){
    
    ## ** Extract quantities from object
    name.endogenous <- endogenous(object)
    n.cluster <- object$data$n

    model.param <- lava::pars(object)
    if(!is.null(param)){
        if(any(names(param) %in% names(model.param) == FALSE)){
            stop("Argument \'param\' have appropriate names: \"",
                 paste(setdiff(names(param),names(model.param)), collapse = "\" \""),
                 "\" \n")
        }
        model.param[names(param)] <- param
    }

    if(is.null(data)){
        data <- object$data$model.frame
    }

    name.param <- names(model.param)

    n.latent <- length(latent(object))
    
    ### ** Compute conditional moments and derivatives
    dMoments <- conditionalMoment(object, data = data, param = model.param,
                                  second.order = df, usefit = TRUE)
    name.meanparam <- names(dMoments$dmu)
    name.varparam <- names(dMoments$dOmega)

    ### ** Compute residuals
    epsilon <- .calcResidualsLVM(data = data, dMoments = dMoments,
                                 n.latent = n.latent,
                                 name.endogenous = name.endogenous)

    ### ** Identify missing values
    if(any(is.na(epsilon))){
        index.Omega <- lapply(1:n.cluster,function(iC){which(!is.na(epsilon[iC,]))})

        ## convert to list format to enable different number of observation per cluster
        dMoments$dmu <- lapply(dMoments$dmu, function(x){
            lapply(1:n.cluster,function(iC){x[iC,index.Omega[[iC]]]})
        })
        dMoments$d2mu <- lapply(dMoments$d2mu, function(x){
            lapply(x, function(y){
                lapply(1:n.cluster,function(iC){y[iC,index.Omega[[iC]]]})
            })
        })
        epsilon <- lapply(1:n.cluster,function(iC){epsilon[iC,index.Omega[[iC]]]})
    }else{
        index.Omega <- NULL
    }

    ### ** Compute residual variance covariance matrix
    Omega <- .calcOmegaLVM(dMoments, n.latent = n.latent)
    
    ## ** param with non-zero third derivative
    type.3deriv <- c("alpha","Gamma","Lambda","B","Psi_var","Sigma_var","Psi_cov","Sigma_cov")
    index.keep <- intersect(which(!is.na(dMoments$df.param$lava)),
                            which(dMoments$df.param$detail %in% type.3deriv)
                            )
    
    name.3deriv <- dMoments$df.param[index.keep, "originalLink"]

    ## ** correction
    if(df == FALSE){
        derivative <- "none"
    }else if(numeric.derivative){
        derivative <- "numeric"
    }else{
        derivative <- "analytic"
    }
    
    out <- .sCorrect(object,
                     param = model.param,
                     epsilon = epsilon,
                     Omega = Omega,
                     dmu = dMoments$dmu,
                     dOmega = dMoments$dOmega,
                     d2mu = dMoments$d2mu,
                     d2Omega = dMoments$d2Omega,
                     name.param = name.param,
                     name.meanparam = name.meanparam,
                     name.varparam = name.varparam,
                     name.endogenous = name.endogenous,
                     name.3deriv = name.3deriv,
                     n.cluster = n.cluster,
                     index.Omega = index.Omega,
                     adjust.Omega = adjust.Omega,
                     adjust.n = adjust.n,
                     tol = tol,
                     n.iter = n.iter,
                     score = score,
                     derivative = derivative,
                     ...)

    ## ** export
    out$args <- list(adjust.Omega = adjust.Omega,
                     adjust.n = adjust.n,                     
                     df = df,
                     numeric.derivative = numeric.derivative,
                     tol = tol, n.iter = n.iter)
    return(out)       
}

## * sCorrect.lvmfit2
#' @rdname sCorrect
#' @export
sCorrect.lvmfit2 <- function(object, ...){
    class(object) <- setdiff(class(object),"lvmfit2")
    return(sCorrect(object, ...))    
}
## * .sCorrect
.sCorrect <- function(object, param, epsilon, Omega, dmu, dOmega, d2mu, d2Omega, 
                      name.param, name.meanparam, name.varparam, name.endogenous, name.3deriv,
                      n.cluster, index.Omega,
                      adjust.Omega, adjust.n, tol, n.iter, score, derivative){

    n.param <- length(param)
    if(!is.null(index.Omega)){
        n.endogenous.cluster <- lapply(index.Omega,length)
    }else{
        n.endogenous.cluster <- NULL
    }
    
    ## ** corrected ML estimates
    out  <- adjustEstimate(epsilon = epsilon,
                           Omega = Omega,
                           dmu = dmu,
                           dOmega = dOmega,
                           n.cluster = n.cluster,
                           name.param = name.param,
                           name.meanparam = name.meanparam,
                           name.varparam = name.varparam,
                           name.endogenous = name.endogenous,
                           index.Omega = index.Omega, ## mode2
                           adjust.Omega = adjust.Omega,
                           adjust.n = adjust.n,
                           tol = tol, n.iter = n.iter)    
    out$param <- param

    browser()
    ## ** corrected score
    if(score){
        out$score <- .score2(epsilon = out$epsilon,
                             Omega = out$Omega,
                             OmegaM1 = out$OmegaM1,
                             dmu = dmu,
                             dOmega = dOmega,
                             name.param = name.param,
                             name.meanparam = name.meanparam,
                             name.varparam = name.varparam,
                             index.Omega = index.Omega, ## mode2
                             n.cluster = n.cluster,
                             indiv = TRUE)
    }
    
    ## ** first derivative of the expected information matrix
    if(derivative == "numeric"){
        if(adjust.Omega || adjust.n){
            warning("The numerical derivative of the information matrix is computed ignoring the small sample correction \n")
        }
        
        ### *** direct computation of the variance-covariance matrix
        calcVcov <- function(iParam){ # x <- p.obj
            pp <- param
            pp[names(iParam)] <- iParam

            args <- object$sCorrect$args
            args$df <- FALSE
            args$score <- FALSE
            vcov.param <- do.call(sCorrect,
                                  args = c(list(object, param = pp), args))$vcov.param
            
            return(vcov.param)
        }

        ### *** numerical derivative
        test.package <- try(requireNamespace("numDeriv"), silent = TRUE)
        if(inherits(test.package,"try-error")){
            stop("There is no package \'numDeriv\' \n",
                 "This package is necessary when argument \'numeric.derivative\' is TRUE \n")
        }
        jac.param <- param[name.3deriv]
        res.numDeriv <- numDeriv::jacobian(calcVcov, x = jac.param, method = "Richardson")
        
        out$dVcov.param <- array(res.numDeriv,
                                 dim = c(n.param,n.param,length(name.3deriv)),
                                 dimnames = list(name.param, name.param, name.3deriv))
        
    }else if(derivative == "analytic"){
        if(is.null(index.Omega)){
            dInfo.dtheta <- .d2Information(dmu = dmu,
                                           d2mu = d2mu,
                                           dOmega = dOmega,
                                           d2Omega = d2Omega,
                                           Omega = out$Omega,
                                           OmegaM1 = out$OmegaM1,
                                           n.corrected = out$n.corrected,
                                           name.param  = name.param,
                                           name.3deriv = name.3deriv)
        }else{
            dInfo.dtheta <- .d2InformationIndiv(dmu = dmu,
                                                d2mu = d2mu,
                                                dOmega = dOmega,
                                                d2Omega = d2Omega,
                                                Omega = out$Omega,
                                                OmegaM1 = out$OmegaM1,
                                                n.corrected = out$n.corrected,
                                                n.endogenous.cluster = n.endogenous.cluster,
                                                index.Omega = index.Omega,
                                                name.param  = name.param,
                                                name.3deriv = name.3deriv)
        }

        p3 <- dim(dInfo.dtheta)[3]
        out$dVcov.param <- array(NA, dim = dim(dInfo.dtheta), dimnames = dimnames(dInfo.dtheta))
        for(iP in 1:p3){
            out$dVcov.param[,,iP] <- - out$vcov.param %*% dInfo.dtheta[,,iP] %*% out$vcov.param 
        }
    }
        
    ## ** export
    return(out)
}

## * sCorrect<-
#' @rdname sCorrect
#' @export
`sCorrect<-` <-
  function(x, ..., value) UseMethod("sCorrect<-")

## * sCorrect<-.lm
#' @rdname sCorrect
#' @export
`sCorrect<-.lm` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    class(x) <- append("lm2",class(x))

    return(x)
}    
## * sCorrect<-.gls
#' @rdname sCorrect
#' @export
`sCorrect<-.gls` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    class(x) <- append("gls2",class(x))

    return(x)
}    
## * sCorrect<-.lme
#' @rdname sCorrect
#' @export
`sCorrect<-.lme` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    class(x) <- append("lme2",class(x))
    
    return(x)
}    

## * sCorrect<-.lvmfit
#' @rdname sCorrect
#' @export
`sCorrect<-.lvmfit` <- function(x, ..., value){
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    class(x) <- append("lvmfit2",class(x))

    return(x)
}    

## * sCorrect<-.lvmfit2
#' @rdname sCorrect
#' @export
`sCorrect<-.lvmfit2` <- function(x, ..., value){

    class(x) <- setdiff(class(x),"lvmfit2")
    x$sCorrect <- sCorrect(x, ..., adjust.Omega = value, adjust.n = value)
    class(x) <- append("lvmfit2",class(x))
    
    return(x)
}


##----------------------------------------------------------------------
### sCorrect.R ends here







## * .calcResidualsLVM
.calcResidualsLVM <- function(data, dMoments, n.latent, name.endogenous){

    ## ** fitted value
    object.fitted <- dMoments$value$nu.XK
    if(n.latent>0){
        object.fitted <- object.fitted + dMoments$value$alpha.XGamma.iIB %*% dMoments$value$Lambda
    }

    ## ** residuals
    out <- data[, name.endogenous] - object.fitted
    
    ## ** export
    return(as.matrix(out))

}

## * .calcOmegaLVM
.calcOmegaLVM <- function(dMoments, n.latent){
    if(n.latent>0){
        Omega <- dMoments$value$tLambda.tiIB.Psi.iIB %*% dMoments$value$Lambda + dMoments$value$Sigma
    }else{
        Omega <- dMoments$value$Sigma
    }
    return(Omega)
}

### residuals2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (09:05) 
## Version: 
## Last-Updated: jan  3 2018 (15:56) 
##           By: Brice Ozenne
##     Update #: 729
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * documentation - residuals2
#' @title Compute the residuals from a lvmfit objects
#' @description Compute the residuals from a lvmfit objects
#' @name residuals2
#' 
#' @param x a fitted latent variable model.
#' @param p [optional] vector of parameters at which to evaluate the score.
#' @param data [optional] data set.
#'
#' @examples
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' e <- estimate(m,sim(m,1e2))
#' residuals2(e)
#' @export
`residuals2` <-
    function(object, ...) UseMethod("residuals2")

## * residuals2.lm

residuals2.lm <- function(object, 
                          adjust.residuals = TRUE,
                          return.vcov.param = FALSE, ...){

### ** Extract information

    ## *** parameters
    p <- coef(object)
    n <- NROW(object$model)
    name.param <- names(p)
    n.param <- length(name.param)
    
    ## *** formula
    formula.object <- formula(object)
    
    ## *** data    
    X <- model.matrix(object)
    
    name.Y <- all.vars(update(formula.object, ".~1"))
    Y <- model.frame(object)[[name.Y]]

### ** Compute observed residuals
    epsilon <- residuals(object, type = "response")

    ### ** Small sample adjustement
    iXX <- solve(t(X) %*% X)
    if(adjust.residuals){
        ## *** Compute the leverage
        leverage <- rowSums((X %*% iXX) * X)
        ### same as influence(object)$hat
        ## range(leverage-influence(object)$hat)
        ## save as (but faster)
        ## leverage - diag(X %*% iXX %*% t(X))        
        if(return.vcov.param){            
            sigma0 <- mean(epsilon^2)
            sigma2 <- sigma0 + mean(sigma0 * leverage)
            ### same as
            ## sigma0 + sigma0 / n * tr(X %*% iXX %*% t(X))
            ## sigma0 * (n+sum(leverage))/n
            ## sigma_corrected / sigma(object)^2
            vcov.param <- sigma2 * iXX            
            ## vcov.param / (1/(n-sum(leverage))*sum(epsilon * epsilon)) * iXX
            ## vcov.param / vcov(object)
            vcov.sigma2 <- 2*sigma2^2/(n-sum(leverage))
        }
        epsilon <- epsilon/sqrt(leverage)
        
    }else if(return.vcov.param){
        sigma2 <- mean(epsilon^2)
        vcov.param <- sigma2 * iXX        
        vcov.sigma2 <- 2*sigma2^2/n 
    }

    if(return.vcov.param){
        vcov.all <- as.matrix(bdiag(vcov.param, vcov.sigma2))
        attr(epsilon, "vcov.param") <- vcov.all
        attr(epsilon, "sigma2") <- sigma2
    }
       
    return(epsilon)
}

## * residuals2.gls
#' @rdname residuals2
#' @export
residuals2.gls <- function(object, cluster = NULL, p = NULL, data = NULL,
                           adjust.residuals = TRUE, power = 1/2, as.clubSandwich = TRUE,
                           second.order = FALSE,
                           return.vcov.param = FALSE, return.prepareScore2 = FALSE, ...){

    if(object$method!="ML"){
        stop("Only work with Maximum Likelihood estimation \n")
    }
    
    test.var <- !is.null(object$modelStruct$varStruct)
    test.cor <- !is.null(object$modelStruct$corStruct)

    if(test.var){
        validClass.var <- c("varIdent","varFunc")
        if(any(class(object$modelStruct$varStruct) %in% validClass.var == FALSE)){
            stop("can only handle varStruct of class \"varIdent\"\n")
        }
    }
    if(test.cor){
        validClass.cor <- c("corCompSymm","corSymm","corStruct")
        if(any(class(object$modelStruct$corStruct) %in% validClass.cor == FALSE)){
            stop("can only handle corStruct of class \"corCompSymm\" or \"corSymm\"\n")
        }
    }
    
### ** Extract information
    ## *** parameters
    pp <- .coef2(object)
    attr.param <- attributes(pp)
    attributes(pp) <- attr.param["names"]

    name.param <- attr.param$names
    n.param <- length(name.param)
    
    ## *** formula
    formula.object <- .getFormula2(object)
    
    ## *** data    
    if(is.null(data)){
        data <- getData(object)
    }
    X <- model.matrix(formula.object, data)
    X <- X[,attr.param$mean.coef,drop=FALSE] ## drop unused columns (e.g. factor with 0 occurence)    
    attr(X,"assign") <- NULL
    attr(X,"contrasts") <- NULL
    
    name.Y <- all.vars(update(formula.object, ".~1"))
    Y <- data[[name.Y]]

    ## *** group
    resGroup <- .getGroups2(object,
                            cluster = cluster,
                            data = data)

    cluster <- resGroup$cluster
    n.cluster <- resGroup$n.cluster
    endogenous <- resGroup$endogenous
    name.endogenous <- resGroup$name.endogenous
    n.endogenous <- resGroup$n.endogenous
    index.obs <- resGroup$index.obs
   
### ** Prepare
    
    
    ## *** update parameters with user-specified values
    if(!is.null(p)){
        
        if(any(name.param %in% names(p)==FALSE)){
            stop("argument \'p\' is not correctly specified \n",
                 "missing parameters: \"",paste(name.param[name.param %in% names(p) == FALSE], collapse = "\" \""),"\"\n")
        }
        
        p <- p[attr.param$names]

    }else{
        
        p <- pp
        
    }    
    
### ** Compute observed residuals
    epsilon <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                      dimnames = list(NULL, name.endogenous))
    epsilon[index.obs] <- Y - X %*% p[attr.param$mean.coef]
    ## residuals(object)-as.vector(t(epsilon))

### ** Reconstruct variance covariance matrix (residuals)
    if(adjust.residuals || return.vcov.param || return.prepareScore2){
        resVcov <- .getVarCov2(object, param = p, attr.param = attr.param,
                               endogenous = endogenous, name.endogenous = name.endogenous, n.endogenous = n.endogenous,
                               cluster = cluster, n.cluster = n.cluster)
    }
    
    ### ** Compute partial derivatives
    if(adjust.residuals || return.vcov.param || return.prepareScore2){
        OPS2 <- prepareScore2(object, X = X,
                              param = p, attr.param = attr.param,
                              n.cluster = n.cluster, name.endogenous = name.endogenous, n.endogenous = n.endogenous,
                              index.obs = index.obs,
                              second.order = second.order)

    }

    ### ** initialize hat matrix
    ls.hat <- lapply(1:n.cluster, function(iC){
            iN.endo <- length(resVcov$ls.indexOmega[[iC]])
            matrix(0, ncol = iN.endo, nrow = iN.endo,
                   dimnames = list(name.endogenous[resVcov$ls.indexOmega[[iC]]],
                                   name.endogenous[resVcov$ls.indexOmega[[iC]]])
                   )
    })
    
    ### ** compute variance covariance matrix (parameters)
    if(adjust.residuals || return.vcov.param){        
        Info <- .information2(dmu.dtheta = OPS2$dmu.dtheta,
                              dOmega.dtheta = OPS2$dOmega.dtheta,
                              Omega = resVcov$Omega,
                              ls.indexOmega = resVcov$ls.indexOmega,
                              hat = ls.hat,
                              n.param = n.param,
                              name.param = name.param,
                              n.cluster = n.cluster)
        vcov.param <- try(chol2inv(chol(Info)), silent = TRUE)
        if("try-error" %in% class(vcov.param)){
            stop("Singular information matrix \n")
        }
        rownames(vcov.param) <- rownames(Info)
        colnames(vcov.param) <- colnames(Info)
        ## vcov.param[rownames(vcov(object)),colnames(vcov(object))] - vcov(object)

        ## factor <- (object$dims$N - object$dims$p)/(object$dims$N - object$dims$p * (object$method == "REML"))
        ## vcov.param[rownames(vcov(object)),colnames(vcov(object))] - vcov(object) * factor
    }
    
    ### ** Normalize residuals    
    if(adjust.residuals){
        resLeverage <- .calcLeverage(dmu.dtheta = OPS2$dmu.dtheta,
                                     dOmega.dtheta = OPS2$dOmega.dtheta,
                                     vcov.param = vcov.param,
                                     Omega = resVcov$Omega, ls.indexOmega = resVcov$ls.indexOmega,
                                     n.cluster = n.cluster, name.endogenous = name.endogenous,
                                     power = power, as.clubSandwich = as.clubSandwich)

        ls.hat <- resLeverage$hat
        epsilon <- do.call(rbind,lapply(1:n.cluster, function(iG){ # iG <- 1
            as.double(resLeverage$powerIH[[iG]] %*% epsilon[iG,])
        }))
        colnames(epsilon) <- name.endogenous
        if(as.clubSandwich<2){
            vcov.param <- resLeverage$vcov.param
            resVcov$Omega <- resLeverage$Omega
        }
     
    }
    
    ### ** export
    if(return.vcov.param){
        attr(epsilon, "vcov.param") <- vcov.param 
    }

    if(return.prepareScore2){
        OPS2$name.param <- name.param
        OPS2$n.param <- n.param
        OPS2$n.cluster <- n.cluster
        OPS2$Omega <- resVcov$Omega
        OPS2$hat <- ls.hat
        OPS2$ls.indexOmega <- resVcov$ls.indexOmega
        attr(epsilon, "prepareScore2") <- OPS2
    }
    return(epsilon)

}

## * residuals2.lme
#' @rdname residuals2
#' @export
residuals2.lme <- residuals2.gls

## * residuals2.lvmfit
#' @rdname residuals2
#' @export
residuals2.lvmfit <- function(object, p = NULL, data = NULL,
                              adjust.residuals = TRUE, power = 1/2, as.clubSandwich = TRUE,
                              second.order = FALSE,
                              return.vcov.param = FALSE, return.prepareScore2 = FALSE,
                              ...){

### ** normalize arguments
    ## test
    if(!identical(class(object),"lvmfit")){
        wrongClass <- paste(setdiff(class(object),"lvmfit"), collapse = " ")
        stop("score2 is not available for ",wrongClass," objects \n")
    }
    if(!is.null(object$model0$attributes$type)){
        stop("score2 is only available for latent variable models involving gaussian variables \n")
    }

    ## param
    if(is.null(p)){
        null.p <- TRUE
        p <- pars(object)
        name.param <- names(p)
    }else{
        null.p <- FALSE
        name.param <- names(pars(object))
        p <- p[name.param]
    }
    n.param <- length(p)

    ## data
    if(is.null(data)){
        data <- model.frame(object)
    }
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }
    n.cluster <- object$data$n
    name.endogenous <- endogenous(object)
    name.latent <- latent(object)
    n.latent <- length(name.latent)
    
    ### ** Reconstruct nu.KX and alpha.XGamma.iIB.Lambda
    if(is.null(object$prepareScore2$update) || object$prepareScore2$update == TRUE){
        OPS2 <- prepareScore2(object, data = data, p = p,
                              name.endogenous = name.endogenous,
                              name.latent = name.latent,
                              second.order = second.order)
    }else{
        OPS2 <- object$prepareScore2
    }

### ** Compute predicted value
    object.fitted <- OPS2$skeleton$value$nu.XK
    if(n.latent>0){
        object.fitted <- object.fitted + OPS2$dtheta$alpha.XGamma.iIB %*% OPS2$skeleton$value$Lambda
    }

### ** Compute residuals
    epsilon <- data[, colnames(object.fitted)] - object.fitted

### ** Compute variance covariance matrix (residuals)
    if(n.latent>0){
        Omega <- OPS2$dtheta$tLambda.tiIB.Psi.iIB %*% OPS2$skeleton$value$Lambda + OPS2$skeleton$value$Sigma
    }else{
        Omega <- OPS2$skeleton$value$Sigma
    }
    ## range(Omega - moments(object, p = p, conditional=TRUE, data = data)$C)

    ### ** Initialize hat matrix
    ls.hat <- lapply(1:n.cluster, function(iC){
        matrix(0, ncol = NCOL(Omega), nrow = NROW(Omega), dimnames = dimnames(Omega))
    })
    
    ### ** Compute variance covariance matrix (parameters)
    if(null.p){
        vcov.param <- vcov(object)
        attr(vcov.param, "det") <- NULL
        attr(vcov.param, "pseudo") <- NULL
        attr(vcov.param, "minSV") <- NULL
    }else{
        Info <- .information2(dmu.dtheta = OPS2$dtheta$dmu.dtheta,
                              dOmega.dtheta = OPS2$dtheta$dOmega.dtheta,
                              Omega = Omega, ls.indexOmega = NULL,
                              hat = ls.hat,
                              n.param = n.param, name.param = name.param, n.cluster = n.cluster)

        vcov.param <- solve(Info)
        rownames(vcov.param) <- rownames(Info)
        colnames(vcov.param) <- colnames(Info)
    }        
    ## round(vcov.param[rownames(vcov(object)),colnames(vcov(object))] - vcov(object),10)

    ### ** Normalize residuals
    if(adjust.residuals){
        resLeverage <- .calcLeverage(dmu.dtheta = OPS2$dtheta$dmu.dtheta,
                                     dOmega.dtheta = OPS2$dtheta$dOmega.dtheta,
                                     vcov.param = vcov.param,
                                     Omega = Omega,
                                     ls.indexOmega = NULL,
                                     n.cluster = n.cluster, name.endogenous = name.endogenous,
                                     power = power, as.clubSandwich = as.clubSandwich)

        epsilon <- do.call(rbind,lapply(1:n.cluster, function(iG){ # iG <- 1
            as.double(resLeverage$powerIH[[iG]] %*% epsilon[iG,])
        }))
        colnames(epsilon) <- name.endogenous
        
        if(as.clubSandwich<2){
            vcov.param <- resLeverage$vcov.param
            Omega <- resLeverage$Omega
        }
        ls.hat <- resLeverage$hat        
        
    }

### ** Export
    if(return.vcov.param){
        attr(epsilon, "vcov.param") <- vcov.param 
    }
    if(return.prepareScore2){
        OPS2$name.param <- name.param
        OPS2$n.param <- n.param
        OPS2$n.cluster <- n.cluster
        OPS2$Omega <- Omega
        OPS2$hat <- ls.hat
        attr(epsilon, "prepareScore2") <- OPS2
    }
    return(epsilon)
}

## * .calcLeverage
.calcLeverage <- function(dmu.dtheta, dOmega.dtheta, vcov.param, 
                          Omega, ls.indexOmega,
                          n.cluster, name.endogenous,
                          power, as.clubSandwich){

    clusterSpecific <- !is.null(ls.indexOmega)
    if(clusterSpecific==FALSE){
        Omega.tempo <- Omega
        Omega_chol.tempo <- chol(Omega)
        iOmega.tempo <- chol2inv(Omega_chol.tempo)
        name.endogenous.tempo <- name.endogenous
        n.endogenous.tempo <- length(name.endogenous.tempo)
    }
    name.param <- colnames(vcov.param)
    n.param <- length(name.param)
    name.meanparam <- names(dmu.dtheta)
    
    
### ** Compute leverage
    ls.biasOmega <- list()
    ls.H <- list()
    ls.powerIH <- list()
    
    for(iG in 1:n.cluster){ ## iG <- 1
        
        ## *** Prepare
        if(clusterSpecific){
            Omega.tempo <- Omega[ls.indexOmega[[iG]],ls.indexOmega[[iG]],drop=FALSE]
            Omega_chol.tempo <- chol(Omega.tempo)
            iOmega.tempo <- chol2inv(Omega_chol.tempo)
            
            name.endogenous.tempo <- name.endogenous[ls.indexOmega[[iG]]]
            n.endogenous.tempo <- length(ls.indexOmega[[iG]])
        }
        Id.tempo <- diag(1,nrow=n.endogenous.tempo,ncol=n.endogenous.tempo)
        
        ## *** Combine derivatives into a matrix
        dmu.dtheta.tempo <- matrix(0, nrow = n.endogenous.tempo, ncol = n.param,
                                   dimnames = list(name.endogenous.tempo, name.param))
        for(iP in name.meanparam){
            dmu.dtheta.tempo[,iP] <- dmu.dtheta[[iP]][iG,name.endogenous.tempo]
        }

        ## *** Compute IH
        ls.biasOmega[[iG]] <- dmu.dtheta.tempo %*% vcov.param %*% t(dmu.dtheta.tempo)
        ls.H[[iG]] <- ls.biasOmega[[iG]] %*% iOmega.tempo
        IH <- Id.tempo - ls.H[[iG]]

        ## correction
        if(power == 1){
            ls.powerIH[[iG]] <- solve(IH)
        }else{
            if(as.clubSandwich){
                M.tempo <- Omega_chol.tempo %*% IH %*% Omega.tempo %*% t(Omega_chol.tempo)
                ls.powerIH[[iG]] <- as.matrix(t(Omega_chol.tempo) %*% matrixPower(M.tempo, power = -1/2, symmetric = TRUE) %*% Omega_chol.tempo)
                ## crossprod(iIH) %*% IH
            }else{
                IH_sym <- iOmega.tempo %*% IH
                iIH_sym <- matrixPower(IH_sym, power = -power, symmetric = TRUE)
                
                ls.powerIH[[iG]] <- chol(iOmega.tempo) %*% iIH_sym
                ## crossprod(iIH_sym) %*% IH_sym
            }
                
        }
    }

    ### ** correct vcov
    Omega <- Omega + Reduce("+",ls.biasOmega)/n.cluster
    
    Info <- .information2(dmu.dtheta = dmu.dtheta,
                          dOmega.dtheta = dOmega.dtheta,
                          Omega = Omega, ls.indexOmega = ls.indexOmega,
                          hat = ls.H,
                          n.param = n.param, name.param = name.param,
                          n.cluster = n.cluster)
    vcov.param <- chol2inv(chol(Info))
    rownames(vcov.param) <- rownames(Info)
    colnames(vcov.param) <- colnames(Info)
    
### ** export
    out <- list(powerIH = ls.powerIH,
                hat = ls.H,
                vcov.param = vcov.param,
                Omega = Omega)
    
    return(out)
}


##----------------------------------------------------------------------
### residuals2.R ends here


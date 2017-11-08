### residuals2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (09:05) 
## Version: 
## Last-Updated: nov  8 2017 (21:47) 
##           By: Brice Ozenne
##     Update #: 235
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

## * residuals2.lvmfit
#' @rdname residuals2
#' @export
residuals2.lvmfit <- function(object, p = NULL, data = NULL,
                              adjust.residuals = TRUE, power = 1/2, as.clubSandwich = TRUE,
                              ls.value = NULL,
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
    name.param.lava <- names(pars(object))
    if(is.null(p)){
        p <- pars(object)        
    }else{
        p <- p[name.param.lava]
    }

    ## data
    if(is.null(data)){
        data <- model.frame(object)
    }
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }
    n <- object$data$n
    name.endogenous <- endogenous(object)
    n.endogenous <- length(name.endogenous)
    name.latent <- latent(object)
    n.latent <- length(name.latent)
   
### ** Reconstruct nu.KX and alpha.XGamma.iIB.Lambda
    if(is.null(ls.value)){   
        if(is.null(object$prepareScore2)){
            ls.value <- skeleton(object, as.lava = TRUE,
                                 n.latent = n.latent, name.latent = name.latent,
                                 n.endogenous = n.endogenous, name.endogenous = name.endogenous,
                                 update.value = TRUE, p = p, data = data)$value
        }else{
            ls.value <- .updateValueWithSkeleton(param = p, data = data,
                                                 n.latent = n.latent, name.latent = name.latent,
                                                 n = n, n.endogenous = n.endogenous, name.endogenous = name.endogenous,
                                                 skeleton = object$prepareScore2$skeleton,
                                                 value = object$prepareScore2$value)
        }
        
        ls.value$alpha.XGamma.iIB.Lambda <- ls.value$alpha.XGamma %*% solve(diag(1,n.latent,n.latent)-ls.value$B) %*% ls.value$Lambda
    }
    
### ** Compute predicted value
    object.fitted <- ls.value$nu.XK
    if(n.latent>0){
        object.fitted <- object.fitted + ls.value$alpha.XGamma.iIB.Lambda
    }

### ** Compute residuals
    epsilon <- data[, colnames(object.fitted)] - object.fitted

### ** Normalize residuals
    if(adjust.residuals){
        ## *** gather derivatives
        dmu.dtheta <- sapply(colnames(vcov.param), function(iP){
            if(is.null(dmu.dtheta[[iP]])){
                return(rep(0, times = n * n.endogenous))
            }else{
                return(as.vector(t(dmu.dtheta[[iP]])))
            }
        })
        colnames(dmu.dtheta) <- colnames(vcov.param)

        ## *** adjust residuals
        ls.leverage <- .calcLeverage(n.group = n, table.group = rep(n.endogenous, n), ls.indexGroup = NULL,
                                     dmu.dtheta = dmu.dtheta, vcov.param = vcov.param,
                                     Omega = Omega, iOmega = iOmega, Omega_chol = Omega_chol,
                                     power = power, as.clubSandwich = as.clubSandwich)

        epsilon <- do.call(rbind,lapply(1:n, function(iG){ # iG <- 1
            as.double(ls.leverage[[iG]] %*% epsilon[iG,])
        }))
    }   
    
### ** Export
    return(epsilon)
}

## * residuals2.gls
#' @rdname residuals2
#' @export
residuals2.gls <- function(object, cluster = NULL, p = NULL, data = NULL,
                           adjust.residuals = TRUE, power = 1/2, as.clubSandwich = TRUE,
                           return.vcov.param = FALSE, return.score = FALSE, ...){

    test.var <- !is.null(object$modelStruct$varStruct)
    test.cor <- !is.null(object$modelStruct$corStruct)
    if(is.null(data)){
        data <- getData(object)
    }
        
### ** Extract information
    ## *** data    
    formula.object <- evalInParentEnv(object$call$model)
    
    X <- model.matrix(formula.object, data)
    attr(X,"assign") <- NULL
    attr(X,"contrasts") <- NULL
    
    name.Y <- all.vars(update(formula.object, ".~1"))
    Y <- data[[name.Y]]

    ## *** mean parameters
    mean.coef <- coef(object)

    ## *** variance parameters
    if(test.var){
        var.coef <- c(sigma2 = sigma(object),coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))^2
    }else{
        var.coef <- c(sigma2 = sigma(object)^2)
    }

    ## *** covariance parameters
    if(test.cor){
        cor.coef <- coef(object$modelStruct$corStruct, unconstrained = FALSE)
        names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
    }else{
        cor.coef <- NULL
    }

    ## *** update with user-specified values
    if(!is.null(p)){
        name.meancoef <- names(mean.coef) 
        name.corcoef <- names(cor.coef)
        name.varcoef <- names(var.coef)
        
        if(any(c(name.meancoef,name.corcoef,name.varcoef) %in% names(p)==FALSE)){
            name.all <- c(name.meancoef,name.corcoef,name.varcoef)
            stop("argument \'p\' is not correctly specified \n",
                 "missing parameters: \"",paste(name.all[name.all %in% names(p) == FALSE], collapse = "\" \""),"\"\n")
        }

        mean.coef <- p[name.meancoef]
        cor.coef <- p[name.corcoef]
        var.coef <- p[name.varcoef]
    }
    name.param <- c(names(mean.coef),names(cor.coef),names(var.coef))
    n.param <- length(name.param)

    ## *** cluster
    if(test.cor){
        cluster <- as.numeric(object$groups)
        if(any(class(object$modelStruct$corStruct) %in% c("corCompSymm","corSymm","corStruct") == FALSE)){
            stop("can only handle corStruct of class \"corCompSymm\" or \"corSymm\"\n")
        }
    }else{
        if(length(cluster) == 1 && is.character(cluster)){
            cluster <- as.numeric(as.factor(data[[cluster]]))
        }else{
            if(length(cluster)!=NROW(data)){
                stop("length of cluster and data do not match \n")
            }
            cluster <- as.numeric(as.factor(cluster))
        }
    }
    if(test.var){
        if(any(class(object$modelStruct$varStruct) %in% c("varIdent","varFunc") == FALSE)){
            stop("can only handle varStruct of class \"varIdent\"\n")
        }
    }
    
### ** Prepare

    ## cluster
    n.cluster <- length(unique(cluster))
    
    ## endogenous
    if(test.var){
        name.endogenous <- attr(object$modelStruct$varStruct,"groupName")
        vec.rep <- attr(object$modelStruct$varStruct,"groups")        
    }else{
        table.cluster <- table(cluster)
        name.endogenous <- 1:max(table.cluster)
        vec.rep <- unlist(lapply(table.cluster,function(iG){1:iG})) ## not good in presence of missing values
    }
    n.endogenous <- length(name.endogenous)
    
    ## convert observations from the vector format to the matrix format
    index.obs <- cluster+(match(vec.rep,name.endogenous)-1)*n.cluster

### ** Compute observed residuals
    epsilon <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                      dimnames = list(NULL, name.endogenous))    
    epsilon[index.obs] <- Y - X %*% mean.coef
    ## residuals(object)-as.vector(t(epsilon))

### ** Compute score
    if(adjust.residuals || return.score){

        ## *** Reconstruct variance covariance matrix (residuals)
        if(test.cor){
            Omega <- lapply(1:n.cluster,function(iC){ # iC <- 1
                M <- unclass(getVarCov(object, individual = iC))
                colnames(M) <- name.endogenous[vec.rep[cluster==iC]]
                rownames(M) <- name.endogenous[vec.rep[cluster==iC]]
                return(M)
            })
        }else{
            if(test.var){
                sigma2.base <- (sigma(object)*coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE))^2
            }else{
                sigma2.base <- rep(sigma(object)^2,n.endogenous)
            }
            Omega <- tapply(vec.rep, cluster, function(iRep){
                M <- diag(sigma2.base[iRep], nrow = length(iRep), ncol = length(iRep))
                colnames(M) <- name.endogenous[vec.rep[cluster==iC]]
                rownames(M) <- name.endogenous[vec.rep[cluster==iC]]
                return(list(M))

            })
        }

        ## *** Compute score
        res.prepare <- prepareScore2(object, X = X, Omega = Omega,
                                     var.coef = var.coef, cor.coef = cor.coef,
                                     n.cluster = n.cluster, name.endogenous = name.endogenous, n.endogenous = n.endogenous,
                                     index.obs = index.obs)
    }

### ** compute variance covariance matrix (parameters)
    if(return.vcov.param || adjust.residuals){
        Info <- .information2(dmu.dtheta = res.prepare$dmu.dtheta,
                              dOmega.dtheta = res.prepare$dOmega.dtheta,
                              iOmega = res.prepare$iOmega,
                              n.param = n.param,
                              name.param = name.param,
                              n.cluster = n.cluster)
        vcov.param <- solve(Info)        
        rownames(vcov.param) <- rownames(Info)
        colnames(vcov.param) <- colnames(Info)        
        ## vcov.param[rownames(vcov(object)),colnames(vcov(object))] / vcov(object) 
    }else{
        vcov.param <- NULL
    }

### ** Normalize residuals    
    if(adjust.residuals){
        ls.leverage <- .calcLeverage(dmu.dtheta = res.prepare$dmu.dtheta, vcov.param = vcov.param,
                                     Omega = Omega, iOmega = res.prepare$iOmega, Omega_chol = res.prepare$Omega_chol,
                                     n.cluster = n.cluster, 
                                     power = power, as.clubSandwich = as.clubSandwich)

        epsilon <- do.call(rbind,lapply(1:n.cluster, function(iG){ # iG <- 1
            as.double(ls.leverage[[iG]] %*% epsilon[iG,])
        }))
        colnames(epsilon) <- name.endogenous
    }   
    
### ** export
    if(return.vcov.param){
        attr(epsilon, "vcov.param") <- vcov.param 
    }
    if(return.score){
        res.prepare$name.param <- name.param
        res.prepare$n.param <- n.param
        res.prepare$n.cluster <- n.cluster
        attr(epsilon, "score") <- res.prepare
    }
    
    return(epsilon)

}

## * residuals2.lme
#' @rdname residuals2
#' @export
residuals2.lme <- function(object, cluster = NULL, p = NULL, data = NULL,
                           adjust.residuals = TRUE, power = 1/2, as.clubSandwich = TRUE,
                           return.vcov.param = FALSE, return.score = FALSE, ...){

    test.var <- !is.null(object$modelStruct$varStruct)
    test.cor <- !is.null(object$modelStruct$corStruct)
    if(is.null(data)){
        data <- getData(object)
    }

### ** Extract information
    ## *** data    
    formula.object <- evalInParentEnv(object$call$fixed)
    
    X <- model.matrix(formula.object, data)
    attr(X,"assign") <- NULL
    attr(X,"contrasts") <- NULL
    
    name.Y <- all.vars(update(formula.object, ".~1"))
    Y <- data[[name.Y]]

    ## *** mean parameters
    mean.coef <- fixef(object)

    ## *** variance parameters
    if(test.var){
        var.coef <- c(sigma2 = sigma(object),coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))^2
    }else{
        var.coef <- c(sigma2 = sigma(object)^2)
    }

    ## *** covariance parameters
    cor.coef <- as.double(getVarCov(object))
    names(cor.coef) <- paste0("corCoef",1:length(cor.coef))

    ## *** update with user-specified values
    if(!is.null(p)){
        name.meancoef <- names(mean.coef) 
        name.corcoef <- names(cor.coef)
        name.varcoef <- names(var.coef)
        
        if(any(c(name.meancoef,name.corcoef,name.varcoef) %in% names(p)==FALSE)){
            name.all <- c(name.meancoef,name.corcoef,name.varcoef)
            stop("argument \'p\' is not correctly specified \n",
                 "missing parameters: \"",paste(name.all[name.all %in% names(p) == FALSE], collapse = "\" \""),"\"\n")
        }

        mean.coef <- p[name.meancoef]
        cor.coef <- p[name.corcoef]
        var.coef <- p[name.varcoef]
    }
    name.param <- c(names(mean.coef),names(cor.coef),names(var.coef))
    n.param <- length(name.param)

    ## *** cluster
    if(test.cor){
        stop("cannot handle lme objects when corStruct is not null \n")
    }
    if(length(getVarCov(object))>1){
        stop("cannot handle lme objects with more than one random effect \n")
    }
    cluster <- as.numeric(object$groups[,1])
    
    if(test.var){
        if(any(class(object$modelStruct$varStruct) %in% c("varIdent","varFunc") == FALSE)){
            stop("can only handle varStruct of class \"varIdent\"\n")
        }
    }
    
### ** Prepare

    ## cluster
    n.cluster <- length(unique(cluster))
    
    ## endogenous
    if(test.var){
        name.endogenous <- attr(object$modelStruct$varStruct,"groupName")
        vec.rep <- attr(object$modelStruct$varStruct,"groups")        
    }else{
        table.cluster <- table(cluster)
        name.endogenous <- 1:max(table.cluster)
        vec.rep <- unlist(lapply(table.cluster,function(iG){1:iG})) ## not good in presence of missing values
    }
    n.endogenous <- length(name.endogenous)
    
    ## convert observations from the vector format to the matrix format
    index.obs <- cluster+(match(vec.rep,name.endogenous)-1)*n.cluster

### ** Compute observed residuals
    epsilon <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                      dimnames = list(NULL, name.endogenous))    
    epsilon[index.obs] <- Y - X %*% mean.coef
    ## residuals(object)-as.vector(t(epsilon))

### ** Compute score
    if(adjust.residuals || return.score){

        ## *** Reconstruct variance covariance matrix (residuals)
        Omega <- unclass(getVarCov(object, individual = 1:n.cluster, type = "marginal"))
        attr(Omega,"group.levels") <- NULL

        ## *** Compute score
        res.prepare <- prepareScore2(object, X = X, Omega = Omega,
                                     var.coef = var.coef, cor.coef = cor.coef,
                                     n.cluster = n.cluster, name.endogenous = name.endogenous, n.endogenous = n.endogenous,
                                     index.obs = index.obs)
    }

### ** compute variance covariance matrix (parameters)
    if(return.vcov.param || adjust.residuals){
        Info <- .information2(dmu.dtheta = res.prepare$dmu.dtheta,
                              dOmega.dtheta = res.prepare$dOmega.dtheta,
                              iOmega = res.prepare$iOmega,
                              n.param = n.param,
                              name.param = name.param,
                              n.cluster = n.cluster)
        vcov.param <- solve(Info)        
        rownames(vcov.param) <- rownames(Info)
        colnames(vcov.param) <- colnames(Info)        
        ## vcov.param[rownames(vcov(object)),colnames(vcov(object))] / vcov(object) 
    }else{
        vcov.param <- NULL
    }

### ** Normalize residuals    
    if(adjust.residuals){
        ls.leverage <- .calcLeverage(dmu.dtheta = res.prepare$dmu.dtheta, vcov.param = vcov.param,
                                     Omega = Omega, iOmega = res.prepare$iOmega, Omega_chol = res.prepare$Omega_chol,
                                     n.cluster = n.cluster, 
                                     power = power, as.clubSandwich = as.clubSandwich)

        epsilon <- do.call(rbind,lapply(1:n.cluster, function(iG){ # iG <- 1
            as.double(ls.leverage[[iG]] %*% epsilon[iG,])
        }))
        colnames(epsilon) <- name.endogenous
    }   
    
### ** export
    if(return.vcov.param){
        attr(epsilon, "vcov.param") <- vcov.param 
    }
    if(return.score){
        res.prepare$name.param <- name.param
        res.prepare$n.param <- n.param
        res.prepare$n.cluster <- n.cluster
        attr(epsilon, "score") <- res.prepare
    }
    
    return(epsilon)

}

## * .calcLeverage
.calcLeverage <- function(dmu.dtheta, vcov.param, 
                          Omega, iOmega, Omega_chol, n.cluster, 
                          power, as.clubSandwich){

    missing <- is.list(iOmega)
    name.param <- colnames(vcov.param)
    n.param <- length(name.param)
    name.meanparam <- names(dmu.dtheta)
    
### ** Compute leverage 
    ls.leverage <- lapply(1:n.cluster, function(iG){ # iG <- 1

        ## *** Prepare
        if(missing){
            iOmega.tempo <- iOmega[[iG]]
            if(power != 1 && as.clubSandwich){            
                Omega_chol.tempo <- Omega_chol[[iG]]
                Omega.tempo <- Omega[[iG]]
            }
        }else{
            iOmega.tempo <- iOmega
            if(power != 1 && as.clubSandwich){            
                Omega_chol.tempo <- Omega_chol
                Omega.tempo <- Omega
            }
        }
        name.endogenous.tempo <- colnames(iOmega.tempo)
        n.endogenous.tempo <- length(name.endogenous.tempo)
        Id.tempo <- diag(1,nrow=n.endogenous.tempo,ncol=n.endogenous.tempo)
        
        ## *** Combine derivatives into a matrix
        dmu.dtheta.tempo <- matrix(0, nrow = n.endogenous.tempo, ncol = n.param,
                                   dimnames = list(name.endogenous.tempo, name.param))
        for(iP in name.meanparam){
            dmu.dtheta.tempo[,iP] <- dmu.dtheta[[iP]][iG,name.endogenous.tempo]
        }
        
        ## *** Compute IH
        IH <- Id.tempo - dmu.dtheta.tempo %*% vcov.param %*% t(dmu.dtheta.tempo) %*% iOmega.tempo

        ## correction
        if(power == 1){
            iIH <- solve(IH)
        }else{
            if(as.clubSandwich){
                M.tempo <- Omega_chol.tempo %*% IH %*% Omega.tempo %*% t(Omega_chol.tempo)
                iIH <- as.matrix(t(Omega_chol.tempo) %*% matrixPower(M.tempo, power = -1/2, symmetric = TRUE) %*% Omega_chol.tempo)
                ## crossprod(iIH) %*% IH
            }else{
                IH_sym <- iOmega.tempo %*% IH
                iIH_sym <- matrixPower(IH_sym, power = -power, symmetric = TRUE)
                
                iIH <- chol(iOmega.tempo) %*% iIH_sym
                ## crossprod(iIH_sym) %*% IH_sym
            }
                
        }
 
        return(iIH)
    })

### ** export
    return(ls.leverage)
}


##----------------------------------------------------------------------
### residuals2.R ends here

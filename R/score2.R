### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: nov  8 2017 (21:36) 
##           By: Brice Ozenne
##     Update #: 1925
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - score2
#' @title Compute the score directly from the gaussian density
#' @description Compute the score directly from the gaussian density
#' @name score2
#' 
#' @param x a linear model or a latent variable model.
#' @param p [optional] vector of parameters at which to evaluate the score.
#' @param data [optional] data set.
#' @param indiv Should the score relative to each observation be exported? Otherwise the total score (i.e. sum over all observations) will be exported.
#' @param adjust.residuals Should the leverage-adjusted residuals be used to compute the score? Otherwise the raw residuals will be used.
#' @param power the exponent used for computing the leverage-adjusted residuals. See the documentation of the \code{\link{iid2}} function for more details.
#' @param as.clubSandwich The method use to apply the \code{power} argument.
#' @param ... not used.
#'
#' @return A matrix.
#' 
#' @details The log-likelihood of a lvm can written:
#' \deqn{
#'   l(\theta|Y,X) = \sum_{i=1}^{n} - \frac{p}{2} log(2\pi) - \frac{1}{2} log|\Omega(\theta))| - \frac{1}{2} (Y_i-\mu_i(\theta)) \Omega^{-1} (Y_i-\mu_i(\theta))^t
#' }
#' So the score is:
#' \deqn{
#'   s(\theta|Y,X) = \sum_{i=1}^{n}
#' - \frac{1}{2} tr(\Omega(\theta)^{-1} \frac{\partial \Omega(\theta)}{\partial \theta})
#' +  \frac{\partial \mu_i(\theta)}{\partial \theta} \Omega^{-1} (Y_i-\mu_i(\theta))^t 
#' + \frac{1}{2} (Y_i-\mu_i(\theta)) \Omega^{-1} \frac{\partial \Omega(\theta)}{\partial \theta} \Omega^{-1} (Y_i-\mu_i(\theta))^t
#' }
#' with:
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial \nu} = 1 }
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial K} = X_i }
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial \alpha} = (1-B)^{-1}\Lambda }
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial \Gamma} = X_i(1-B)^{-1}\Lambda }
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial \lambda} = (\alpha + X_i \Gamma)(1-B)^{-1}\delta_{\lambda \in \Lambda} }
#' \deqn{ \frac{\partial \mu_i(\beta)}{\partial b} = (\alpha + X_i \Gamma)(1-B)^{-1}\delta_{b \in B}(1-B)^{-1}\Lambda }
#' and:
#' \deqn{ \frac{\partial \Omega(\beta)}{\partial \psi} = \Lambda^t (1-B)^{-t} \delta_{\psi \in \Psi} (1-B) \Lambda }
#' \deqn{ \frac{\partial \Omega(\beta)}{\partial \sigma} = \delta_{\sigma \in \Sigma} }
#' \deqn{ \frac{\partial \Omega(\beta)}{\partial \lambda} = \delta_{\lambda \in \Lambda} (1-B)^{-t} \Psi (1-B)^{-1} \Lambda + \Lambda^t (1-B)^{-t} \Psi (1-B)^{-1} \delta_{\lambda \in \Lambda} }
#' \deqn{ \frac{\partial \Omega(\beta)}{\partial b} = \Lambda^t (1-B)^{-t} \delta_{b \in B} (1-B)^{-t} \Psi (1-B)^{-1} \Lambda - \Lambda^t (1-B)^{-t} \Psi (1-B)^{-1} \delta_{b \in B} (1-B)^{-1} \Lambda}
#'
#' @examples
#'
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' e <- estimate(m,sim(m,1e2))
#' score2(e)
#' 
#' @export
`score2` <-
  function(x, ...) UseMethod("score2")

## * score2.gls
#' @rdname score2
#' @export
score2.gls <- function(object, cluster, p = NULL, data = NULL,
                       adjust.residuals = TRUE, power = 1/2,
                       indiv = TRUE, as.clubSandwich = TRUE, return.vcov.param = FALSE, ...){

### ** Compute residuals (and score)
    epsilon <- residuals2(object, cluster = cluster, p = p, data = data,
                          adjust.residuals = adjust.residuals,
                          power = power,
                          as.clubSandwich = as.clubSandwich,
                          return.vcov.param = return.vcov.param,
                          return.score = TRUE)
    
    res.prepare <- attr(epsilon, "score")
    vcov.param <- attr(epsilon, "vcov.param")
    attr(epsilon, "score") <- NULL
    attr(epsilon, "vcov.param") <- NULL
        
### ** Compute score
    out.score <- .score2(dmu.dtheta = res.prepare$dmu.dtheta,
                         dOmega.dtheta = res.prepare$dOmega.dtheta,
                         epsilon = epsilon,
                         iOmega = res.prepare$iOmega,                         
                         indiv = indiv,
                         name.param = res.prepare$name.param,
                         n.param = res.prepare$n.param,
                         n.cluster = res.prepare$n.cluster)     
    
### ** Export
    if(return.vcov.param){
        attr(out.score,"vcov.param") <- vcov.param
    }
    return(out.score)
}

## * score2.lme
#' @rdname score2
#' @export
score2.lme <- score2.gls

## * score2.lvmfit
#' @rdname score2
#' @export
score2.lvmfit <- function(object, p = NULL, data = NULL, 
                          adjust.residuals = TRUE, power = 1/2,
                          indiv = TRUE, as.clubSandwich = TRUE, return.vcov.param = FALSE, ...){

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
        null.p <- TRUE
        p <- pars(object)        
    }else{
        null.p <- FALSE
        p <- p[name.param.lava]
    }

    ## data
    if(is.null(data)){
        data <- model.frame(object)
    }
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }
    n <- NROW(data)
    
### ** Prepare elements
    name.latent <- latent(object)
    n.latent <- length(name.latent)

    name.endogenous <- endogenous(object)
    n.endogenous <- length(name.endogenous)

    if(is.null(object$prepareScore2)){   
        object$prepareScore2 <- prepareScore2(object, data)
    }
    name.param <- names(object$prepareScore2$type)
    n.param <- length(name.param)
    
### ** Reconstruct Sigma, Lambda, B, Psi
    ls.value <- .updateValueWithSkeleton(param = p, data = data,
                                         n.latent = n.latent, name.latent = name.latent,
                                         n = n, n.endogenous = n.endogenous, name.endogenous = name.endogenous,
                                         skeleton = object$prepareScore2$skeleton,
                                         value = object$prepareScore2$value)

### ** Compute partial derivatives regarding the mean and the variance
    if(any(object$prepareScore2$toUpdate)){
        type2update <- object$prepareScore2$type[object$prepareScore2$toUpdate]
            
        iIB <- solve(diag(1,n.latent,n.latent)-ls.value$B)
        alpha.XGamma.iIB <- ls.value$alpha.XGamma %*% iIB
        iIB.Lambda <-  iIB %*% ls.value$Lambda    
        Psi.iIB <- ls.value$Psi %*% iIB
        tLambda.tiIB.Psi.iIB <- t(iIB.Lambda) %*% Psi.iIB
        
        ## *** mean parameters
        type.meanparam <- type2update[type2update %in% c("alpha","Lambda","Gamma","B")]
        n.meanparam <- length(type.meanparam)
        name.meanparam <- names(type.meanparam)

        if(n.meanparam>0){
            for(iP in 1:n.meanparam){ # iP <- 1
                iType <- type.meanparam[iP]
                iName <- name.meanparam[iP]
            
                if(iType == "alpha"){
                    object$prepareScore2$dmu.dtheta[[iName]] <- object$prepareScore2$dmu.dtheta[[iName]] %*% iIB.Lambda
                }else if(iType == "Gamma"){
                    object$prepareScore2$dmu.dtheta[[iName]] <- object$prepareScore2$dmu.dtheta[[iName]] %*% iIB.Lambda 
                }else if(iType == "Lambda"){
                    object$prepareScore2$dmu.dtheta[[iName]] <- alpha.XGamma.iIB %*% object$prepareScore2$dLambda.dtheta[[iName]]
                }else if(iType == "B"){
                    object$prepareScore2$dmu.dtheta[[iName]] <- alpha.XGamma.iIB %*% object$prepareScore2$dB.dtheta[[iName]] %*% iIB.Lambda
                }

                colnames(object$prepareScore2$dmu.dtheta[[iName]]) <- name.endogenous
            }
        }

        ## *** variance-covariance parameters
        type.vcovparam <- type2update[type2update %in% c("Psi_var","Psi_cov","Lambda","B")]
        n.vcovparam <- length(type.vcovparam)
        name.vcovparam <- names(type.vcovparam)

        if(n.vcovparam>0){
            for(iP in 1:n.vcovparam){ # iP <- 1
                iType <- type.vcovparam[iP]
                iName <- name.vcovparam[iP]
        
                if(iType %in% "Psi_var"){
                    object$prepareScore2$dOmega.dtheta[[iName]] <-  t(iIB.Lambda) %*% object$prepareScore2$dPsi.dtheta[[iName]] %*% iIB.Lambda
                }else if(iType %in% "Psi_cov"){
                    object$prepareScore2$dOmega.dtheta[[iName]] <-  t(iIB.Lambda) %*% object$prepareScore2$dPsi.dtheta[[iName]] %*% iIB.Lambda
                }else if(iType == "Lambda"){
                    object$prepareScore2$dOmega.dtheta[[iName]] <- tLambda.tiIB.Psi.iIB %*% object$prepareScore2$dLambda.dtheta[[iName]]
                    object$prepareScore2$dOmega.dtheta[[iName]] <- object$prepareScore2$dOmega.dtheta[[iName]] + t(object$prepareScore2$dOmega.dtheta[[iName]])
                }else if(iType == "B"){
                    object$prepareScore2$dOmega.dtheta[[iName]] <- tLambda.tiIB.Psi.iIB %*% object$prepareScore2$dB.dtheta[[iName]] %*% iIB.Lambda
                    object$prepareScore2$dOmega.dtheta[[iName]] <- object$prepareScore2$dOmega.dtheta[[iName]] + t(object$prepareScore2$dOmega.dtheta[[iName]])
                }

                colnames(object$prepareScore2$dOmega.dtheta[[iName]]) <- name.endogenous
                rownames(object$prepareScore2$dOmega.dtheta[[iName]]) <- name.endogenous
            }
        }
        
    }

### ** Reconstruct variance covariance matrix (residuals)
    if(n.latent>0){
        Omega <- tLambda.tiIB.Psi.iIB %*% ls.value$Lambda + ls.value$Sigma
    }else{
        Omega <- ls.value$Sigma
    }
    Omega_chol <- chol(Omega)
    iOmega <- chol2inv(Omega_chol) ## faster compared to solve
    ## range(Omega - moments(object, p = p, conditional=TRUE, data = data)$C)

### ** Compute variance covariance matrix (parameters)
    if(adjust.residuals || return.vcov.param){
        if(null.p){
            vcov.param <- vcov(object)
            attr(vcov.param, "det") <- NULL
            attr(vcov.param, "pseudo") <- NULL
            attr(vcov.param, "minSV") <- NULL
        }else{
            Info <- .information2(object$prepareScore2$dmu.theta, object$prepareScore2$dOmega.dtheta,
                                  iOmega = iOmega,
                                  n.param = n.param, name.param = name.param)
            vcov.param <- chol2inv(chol(Info))
            rownames(vcov.param) <- rownames(Info)
            colnames(vcov.param) <- colnames(Info)
        }        
    }else{
        vcov.param <- NULL
    }
    ## round(vcov.param[rownames(vcov(object)),colnames(vcov(object))] - vcov(object),10)

### ** Compute residuals
    epsilon <- residuals2(object, p = p, data = data,
                          adjust.residuals = adjust.residuals,
                          power = power,
                          as.clubSandwich = as.clubSandwich,
                          ls.value = list(vcov.param = vcov.param,
                                          Omega = Omega,
                                          Omega_chol = Omega_chol,
                                          nu.XK = ls.value$nu.XK,
                                          alpha.XGamma.iIB.Lambda = alpha.XGamma.iIB %*% ls.value$Lambda))
    
### ** Compute score
    out.score <- .score2(dmu.dtheta = object$prepareScore2$dmu.dtheta,
                         dOmega.dtheta = object$prepareScore2$dOmega.dtheta,
                         epsilon = as.matrix(epsilon),
                         iOmega = iOmega,
                         indiv = indiv,
                         name.param = name.param, n.param = n.param, n = n)        

### ** export
    out.score <- out.score[,name.param.lava,drop=FALSE]
    if(return.vcov.param){
        attr(out.score,"vcov.param") <- vcov.param[name.param.lava,name.param.lava,drop=FALSE]
    }
    return(out.score)

}

## * .information2
.information2 <- function(dmu.dtheta, dOmega.dtheta, iOmega,
                          n.param, name.param, n.cluster){

    missing <- is.list(iOmega)    
                                        
### ** Compute information
    Info <- matrix(NA, nrow = n.param, ncol = n.param, dimnames = list(name.param,name.param))
    
    for(iP1 in 1:n.param){ # iP <- 1
        for(iP2 in iP1:n.param){ # iP <- 1
            iName1 <- name.param[iP1]
            iName2 <- name.param[iP2]

            if(!is.null(dmu.dtheta[[iName1]]) && !is.null(dmu.dtheta[[iName2]])){
                if(missing){
                    ls.term1 <- lapply(1:n.cluster, function(iC){ # iC <- 1
                        dmu.dtheta[[iName1]][iC,,drop=FALSE] %*% iOmega[[iC]] * dmu.dtheta[[iName2]][iC,]
                    })
                    term1 <- sum(unlist(ls.term1), na.rm = TRUE)
                }else{
                    term1 <- sum(dmu.dtheta[[iName1]] %*% iOmega * dmu.dtheta[[iName2]])
                }
            }else{
                term1 <- 0
            }

            if(!is.null(dOmega.dtheta[[iName1]]) && !is.null(dOmega.dtheta[[iName2]])){
                if(missing){
                    ls.term2 <- lapply(1:n.cluster, function(iC){
                        1/2*tr(iOmega[[iC]] %*% dOmega.dtheta[[iName1]][[iC]] %*% iOmega[[iC]] %*% dOmega.dtheta[[iName2]][[iC]])
                    })
                    term2 <- sum(unlist(ls.term2), na.rm = TRUE)
                }else{                    
                    term2 <- n.cluster/2*tr(iOmega %*% dOmega.dtheta[[iName1]] %*% iOmega %*% dOmega.dtheta[[iName2]])
                }
            }else{
                term2 <- 0
            }
            Info[iP1,iP2] <- (term1+term2)
        }
    }
    Info <- symmetrize(Info, update.upper = FALSE)

### ** export
    return(Info)
}

## * .score2
.score2 <- function(dmu.dtheta, dOmega.dtheta, epsilon, iOmega,
                    indiv,
                    name.param, n.param, n.cluster){
  
    missing <- is.list(iOmega)
    if(missing == FALSE){
        epsilon.iOmega <- epsilon %*% iOmega
    }else{
        epsilon.iOmega <- lapply(1:n.cluster, function(iC){
            name.endogenous.tempo <- colnames(iOmega[[iC]])
            iOmega[[iC]] %*% cbind(epsilon[iC,name.endogenous.tempo])
        })
    }

    name.meanparam <- names(dmu.dtheta)
    name.vcovparam <- names(dOmega.dtheta)
    out.score <- matrix(0, nrow = n.cluster, ncol = n.param, dimnames = list(NULL,name.param))
    
### ** Compute score relative to the mean parameters
    for(iP in name.meanparam){ # iP <- 1
        if(missing){
            term1 <- sapply(1:n.cluster, function(iC){ # iC <- 1
                name.endogenous.tempo <- colnames(iOmega[[iC]])
                dmu.dtheta[[iP]][iC,name.endogenous.tempo] %*% epsilon.iOmega[[iC]]
            })                   
            out.score[,iP] <- out.score[,iP] + term1
        }else{
            out.score[,iP] <- out.score[,iP] + rowSums(dmu.dtheta[[iP]] * epsilon.iOmega) 
        }
            
    }
    
### ** Compute score relative to the variance-covariance parameters
    for(iP in name.vcovparam){ # iP <- 1
        if(missing){
            term2 <- - 1/2 * sapply(1:n.cluster, function(iC){
                tr(iOmega[[iC]] %*% dOmega.dtheta[[iP]][[iC]])
                })
            term3 <- 1/2 * sapply(1:n.cluster, function(iC){
                sum(epsilon.iOmega[[iC]] * dOmega.dtheta[[iP]][[iC]] %*% epsilon.iOmega[[iC]])
                })
        }else{
            term2 <- - 1/2 * tr(iOmega %*% dOmega.dtheta[[iP]])            
            term3 <- 1/2 * rowSums(epsilon.iOmega %*% dOmega.dtheta[[iP]] * epsilon.iOmega)           
        }
        out.score[,iP] <- out.score[,iP] + as.double(term2) + term3 
    }
        

### ** export
    if(indiv==FALSE){
        out.score <- colSums(out.score)
    }
    return(out.score)
}




#----------------------------------------------------------------------
### score2.R ends here

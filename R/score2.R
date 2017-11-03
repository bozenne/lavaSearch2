### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: nov  3 2017 (13:38) 
##           By: Brice Ozenne
##     Update #: 1444
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
#' @export
`score2` <-
  function(x, ...) UseMethod("score2")

## * score2.gls
score2.gls <- function(x, cluster, p = NULL,  data = NULL,
                       adjust.residuals = TRUE, power = 1/2,
                       indiv = TRUE, as.clubSandwich = TRUE, return.vcov.param = FALSE, ...){

### ** Normalize arguments
    if(!identical(class(x),"gls")){
        wrongClass <- paste(setdiff(class(x),"gls"), collapse = " ")
        stop("iid2 is not available for ",wrongClass," objects \n")
    }
    test.var <- !is.null(x$modelStruct$varStruct)
    test.cor <- !is.null(x$modelStruct$corStruct)
    if(test.var == FALSE && test.cor == FALSE){
        stop("Either the correlation or weight argument must be specified when fitting the gls model \n")
    }

    if(is.null(data)){
        data <- getData(x)        
    }
    X <- model.matrix(formula(x), data)

    if(test.cor){
        vec.group <- as.numeric(x$groups)
        if(all(diff(vec.group)>=0) == FALSE && all(diff(vec.group<=0)) == FALSE){
            stop("the grouping variable for the random effect must be sorted before fitting the model \n")
        }
    }else{
        if(length(cluster) == 1 && is.character(cluster)){
            vec.group <- as.numeric(as.factor(data[[cluster]]))
        }else{
            if(length(cluster)!=NROW(data)){
                stop("length of cluster and data do not match \n")
            }
            vec.group <- as.numeric(as.factor(cluster))
        }
    }

    name.coef <- names(coef(x))
    n.coef <- length(name.coef)
    if(!is.null(p)){
        p <- p[name.coef]
    }
    
### ** compute residuals
    if(is.null(p)){
        epsilon <- residuals(x, type = "response", level = 0)
    }else{
        epsilon <- x$fitted+x$residuals - model.matrix(formula(x), getData(x)) %*% p
    }

### ** identify groups
    table.group <- table(vec.group)
    Ugroup <- unique(vec.group)
    n.group <- length(Ugroup)

    ls.indexGroup <- lapply(1:n.group, function(iG){
        which(vec.group==iG)
    })

### ** Reconstruct varaince covariance matrix (residuals)
    if(test.cor){
        tau <- unclass(getVarCov(x))
        ls.Omega <- lapply(1:n.group, function(iI){tau})        
    }else{
        sigma2.base <- (sigma(x)*coef(x$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE))^2
        vec.sigma2 <- sigma2.base[match(attr(x$modelStruct$varStruct,"groups"), names(sigma2.base))]
        
        ls.Omega <- lapply(1:n.group, function(iG){ # iG <- 1
            n.tempo <- table.group[iG]        
            sigma.tempo <- vec.sigma2[ls.indexGroup[[iG]]]
            return(diag(sigma.tempo, n.tempo, n.tempo))
        })
    }
    ls.iOmega <- lapply(ls.Omega, solve)

### ** compute variance covariance matrix (parameters)
    ## factor <- (x$dims$N - (x$method == "REML") * x$dims$p)/(x$dims$N-x$dims$p)
    ## vcov.param <- vcov(x) / factor
    ## check:
    vcov.param <- solve(t(X) %*% Matrix::bdiag(ls.iOmega) %*% X)
                   

### ** Compute score
    score0 <- .score2_nlme(name.param = name.coef,
                           dmu.dtheta = X, epsilon = epsilon,
                           vcov.param = vcov.param, ls.Omega = ls.Omega, ls.iOmega = ls.iOmega,
                           ls.indexGroup = ls.indexGroup, n.group = n.group, table.group = table.group,
                           indiv = indiv, adjust.residuals = adjust.residuals, power = power, as.clubSandwich = as.clubSandwich)
    
### ** Export
    if(return.vcov.param){
        attr(score0,"vcov.param") <- vcov.param
    }
    return(score0)
}


## * score2.lme
score2.lme <- function(x, p = NULL, data = NULL,
                       adjust.residuals = TRUE, power = 1/2,
                       indiv = TRUE, as.clubSandwich = TRUE, return.vcov.param = FALSE, ...){

### ** Normalize arguments
    if(!identical(class(x),"lme")){
        wrongClass <- paste(setdiff(class(x),"lme"), collapse = " ")
        stop("iid2 is not available for ",wrongClass," objects \n")
    }
    if(!is.null(x$modelStruct$corStruct)){
        stop("cannot handle lme objects when corStruct is not null \n")
    }
    if(length(getVarCov(x))>1){
        stop("cannot handle lme objects with more than one random effect \n")
    }
    vec.group <- as.numeric(x$groups[,1])
    if(all(diff(vec.group)>=0) == FALSE && all(diff(vec.group<=0)) == FALSE ){
        stop("the grouping variable for the random effect must be sorted before fitting the model \n")
    }
   
    if(is.null(data)){
        data <- getData(x)        
    }
    X <- model.matrix(formula(x), data)

    name.coef <- names(fixef(x))
    n.coef <- length(name.coef)
    if(!is.null(p)){
        p <- p[name.coef]
    }

    
### ** compute residuals
    if(is.null(p)){
        epsilon <- residuals(x, type = "response", level = 0)
    }else{
        epsilon <- x$fitted[,"fixed"]+x$residuals[,"fixed"] - model.matrix(formula(x), getData(x)) %*% p
        # epsilon - residuals(x, type = "response", level = 0)
    }

### ** identify groups
    table.group <- table(vec.group)
    Ugroup <- unique(vec.group)
    n.group <- length(Ugroup)

    ls.indexGroup <- lapply(1:n.group, function(iG){
        which(vec.group==iG)
    })
    
### ** Reconstruct covariance matrix
    tau <- as.double(getVarCov(x))
    sigma2.base <- sigma(x)^2
    if(!is.null(x$modelStruct$varStruct)){
        Usigma2 <- sigma2.base*coef(x$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)^2
        vec.sigma2 <- Usigma2[match(attr(x$modelStruct$varStruct,"groups"), names(Usigma2))]
    }else{
        vec.sigma2 <- rep(sigma2.base, NROW(data))
    }
    
    ls.Omega <- lapply(1:n.group, function(iG){ # iG <- 1
        n.tempo <- table.group[iG]        
        sigma.tempo <- vec.sigma2[ls.indexGroup[[iG]]]
        return(matrix(tau,nrow=n.tempo,ncol=n.tempo) + diag(sigma.tempo, n.tempo, n.tempo))
    })
    ls.iOmega <- lapply(ls.Omega, solve)
    
### ** compute variance covariance matrix (parameters)
    vcov.param <- solve(t(X) %*% Matrix::bdiag(ls.iOmega) %*% X)
    ## check
    ## round(vcov.param - vcov(x),10)
    
    
### ** Compute score
    score0 <- .score2_nlme(name.param = name.coef,
                           dmu.dtheta = X, epsilon = epsilon,
                           vcov.param = vcov.param, ls.Omega = ls.Omega, ls.iOmega = ls.iOmega,
                           ls.indexGroup = ls.indexGroup, n.group = n.group, table.group = table.group,
                           indiv = indiv, adjust.residuals = adjust.residuals, power = power, as.clubSandwich = as.clubSandwich)

### ** Export
    if(return.vcov.param){
        attr(score0,"vcov.param") <- vcov.param
    }
    return(score0)
}

## * score2.lvmfit
score2.lvmfit <- function(x, p = NULL, data = NULL, 
                          adjust.residuals = TRUE, Dmethod = FALSE, power = 1/2,
                          indiv = TRUE, as.clubSandwich = TRUE, return.vcov.param = FALSE, ...){

### ** normalize arguments
    if(!identical(class(x),"lvmfit")){
        wrongClass <- paste(setdiff(class(x),"lvmfit"), collapse = " ")
        stop("score2 is not available for ",wrongClass," objects \n")
    }
    if(!is.null(x$model0$attributes$type)){
        stop("score2 is only available for latent variable models involving gaussian variables \n")
    }
    
    name.param.lava <- names(pars(x))
    if(is.null(p)){
        null.p <- TRUE
        p <- pars(x)        
    }else{
        null.p <- FALSE
        p <- p[name.param.lava]
    }
   
    if(is.null(data)){
        data <- model.frame(x)
    }
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }
    
### ** Prepare elements
    name.latent <- latent(x)
    n.latent <- length(name.latent)
    
    if(is.null(x$prepareScore2)){   
        x.score2 <- prepareScore2(x, data)
    }else{
        x.score2 <- x$prepareScore2
    }
    name.param <- names(x.score2$type)
    n.param <- length(name.param)
    
    

### ** Reconstruct Sigma, Lambda, B, Psi
    x.score2$Sigma[!is.na(x.score2$skeleton$Sigma)] <- p[x.score2$skeleton$Sigma[!is.na(x.score2$skeleton$Sigma)]]
    if(n.latent>0){
        x.score2$alpha[!is.na(x.score2$skeleton$alpha)] <- p[x.score2$skeleton$alpha[!is.na(x.score2$skeleton$alpha)]]
        for(iLatent in 1:n.latent){
            if(length(x.score2$Gamma[[iLatent]])>0){
                x.score2$Gamma[[iLatent]][!is.na(x.score2$skeleton$Gamma[[iLatent]])] <- p[x.score2$skeleton$Gamma[[iLatent]][!is.na(x.score2$skeleton$Gamma[[iLatent]])]]
            }
        }

        x.score2$Lambda[!is.na(x.score2$skeleton$Lambda)] <- p[x.score2$skeleton$Lambda[!is.na(x.score2$skeleton$Lambda)]]
        x.score2$B[!is.na(x.score2$skeleton$B)] <- p[x.score2$skeleton$B[!is.na(x.score2$skeleton$B)]]
        x.score2$Psi[!is.na(x.score2$skeleton$Psi)] <- p[x.score2$skeleton$Psi[!is.na(x.score2$skeleton$Psi)]]
        
        x.score2$iIB <- solve(diag(1,n.latent,n.latent)-x.score2$B)
        x.score2$iIB.Lambda <-  x.score2$iIB %*% x.score2$Lambda
        x.score2$Psi.iIB <- x.score2$Psi %*% x.score2$iIB
        x.score2$tLambda.tiIB.Psi.iIB <- t(x.score2$iIB.Lambda) %*% x.score2$Psi.iIB
        x.score2$alpha.GammaX <- matrix(NA,nrow = n, ncol = n.latent, byrow = TRUE,
                                        dimnames = list(NULL,name.latent))
        for(iLatent in 1:n.latent){

            if(length(x.score2$Gamma[[iLatent]])>0){
                x.score2$alpha.GammaX[,iLatent] <- x.score2$alpha[iLatent] + data[,x.score2$skeleton$XGamma[[iLatent]],drop=FALSE] %*% x.score2$Gamma[[iLatent]]
            }else{
                x.score2$alpha.GammaX[,iLatent] <- x.score2$alpha[iLatent]
            }
        }
        x.score2$alpha.GammaX.iIB <- x.score2$alpha.GammaX %*% x.score2$iIB

    }

    
   
### ** Reconstruct variance covariance matrix (residuals)
    if(n.latent>0){
        Omega <- x.score2$tLambda.tiIB.Psi.iIB %*% x.score2$Lambda + x.score2$Sigma
    }else{
        Omega <- x.score2$Sigma
    }
    ## range(Omega - moments(x, p = p, conditional=TRUE, data = data)$C)
    Omega_chol <- chol(Omega)
    iOmega <- chol2inv(Omega_chol) ## faster compared to solve

### ** Compute variance covariance matrix (parameters)
    if(adjust.residuals || return.vcov.param){
        if( (null.p || all(abs(pars(x)-p)<1e-10)) && Dmethod == FALSE  ){
            vcov.param <- vcov(x)
        }else{
            ## could be replaced by explicit computation by computing the second derivative of the likelihood
            vcov.param <- -lava::Inverse(numDeriv::jacobian(function(iP){score(x, p = iP, indiv = FALSE, ...)},
                                                            p, method = Dmethod))
            colnames(vcov.param) <- name.param
            rownames(vcov.param) <- name.param
        }        
    }
    ## vcov.param - vcov(x)

### ** Compute partial derivatives
    if(any(x.score2$toUpdate)){
        x.score2 <- .calcDtheta(x.score2)        
    }
    
                                     
### ** Residuals
    epsilon <- as.matrix(residuals(x, p = p))
    
    if(adjust.residuals){

        ## *** gather derivatives
        dmu.dtheta <- sapply(colnames(vcov.param), function(iP){
            if(is.null(x.score2$dmu.dtheta[[iP]])){
                return(rep(0, times = n * n.endogenous))
            }else{
                return(as.vector(t(x.score2$dmu.dtheta[[iP]])))
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
    }else{
        ls.leverage <- NULL
    }   
   
### ** Compute score
    n.score <- indiv*n+(1-indiv)*1 ## n if indiv==TRUE and 1 if indiv==FALSE
    out.score <- matrix(NA, nrow = n.score, ncol = n.param, dimnames = list(NULL,name.param))
    fct.sum <- switch(as.character(as.double(indiv)),
                      "1" = "rowSums",
                      "0" = "sum")

    ## tiOmega.epsilon <- t(iOmega %*% t(epsilon))
    tiOmega.epsilon <- epsilon %*% iOmega

    for(iP in 1:n.param){ # iP <- 1
        iScore <- rep(0,n.score)
        iName <- name.param[iP]
        if(!is.null(x.score2$dmu.dtheta[[iName]])){
            iScore <- do.call(fct.sum, args = list(x.score2$dmu.dtheta[[iName]] * tiOmega.epsilon))
        }

        if(!is.null(x.score2$dOmega.dtheta[[iName]])){
            firstTerm <- - (indiv*1+(1-indiv)*n)/2 * tr(iOmega %*% x.score2$dOmega.dtheta[[iName]])            
            secondTerm <- 1/2 * do.call(fct.sum, args = list((tiOmega.epsilon %*% x.score2$dOmega.dtheta[[iName]]) * tiOmega.epsilon))
            iScore <- iScore + as.double(firstTerm) + secondTerm            
        }
        out.score[,iP] <- iScore       
    }

### ** export
    if(return.vcov.param){
        attr(out.score,"vcov.param") <- vcov.param
    }
    return(out.score[,name.param.lava,drop=FALSE])

}

## * .score2_nlme
.score2_nlme <- function(name.param,
                         dmu.dtheta, epsilon, vcov.param, ls.Omega, ls.iOmega,
                         indiv, adjust.residuals, power, as.clubSandwich,
                         ls.indexGroup, n.group, table.group){

    Omega <- Matrix::bdiag(ls.Omega)
    iOmega <- Matrix::bdiag(ls.iOmega)
    Omega_chol <- Matrix::bdiag(lapply(ls.Omega,chol))
    
    ## ** compute leverage-adjusted residuals
    if(adjust.residuals){
        ls.leverage <- .calcLeverage(n.group = n.group, table.group = table.group, ls.indexGroup = ls.indexGroup,
                                     dmu.dtheta = dmu.dtheta, vcov.param = vcov.param,
                                     Omega = Omega, iOmega = iOmega, Omega_chol = Omega_chol,
                                     power = power, as.clubSandwich = as.clubSandwich)

        epsilon <- do.call(rbind,lapply(1:n.group, function(iG){ # iG <- 1
            ls.leverage[[iG]] %*% epsilon[ls.indexGroup[[iG]]]
        }))
        
    }else{

        ls.leverage <- NULL
        
    }
    
    ## ** Compute score
    score0 <- t(sapply(1:n.group, function(iG){ # iG <- 1
        index.tempo <- ls.indexGroup[[iG]]
        score.tempo <- t(dmu.dtheta[index.tempo,,drop=FALSE]) %*% ls.iOmega[[index.tempo]] %*% epsilon[index.tempo]
        return(as.double(score.tempo))
    }))
   
    ## ** export
    if(indiv==FALSE){
        score0 <- rbind(colSums(score0))
    }
    colnames(score0) <- name.param

    return(score0)
   
}


## * .calcDetheta
.calcDtheta <- function(x){

    type <- x$type[x$toUpdate]
    param <- x$param[x$toUpdate]
    n.param <- length(param)
    
### ** partial derivative

    for(iP in 1:n.param){ # iP <- 1
        iType <- type[iP]
        iName <- param[iP]

        cat(iName,": ", iType,"\n")

        ## *** linear predictor
        if(iType == "alpha"){
            x$dmu.dtheta[[iName]] <- x$dmu.dtheta[[iName]] %*% x$iIB.Lambda            
        }else if(iType == "Gamma"){
            x$dmu.dtheta[[iName]] <- x$dmu.dtheta[[iName]] %*% x$iIB.Lambda 
        }else if(iType == "Lambda"){
            x$dmu.dtheta[[iName]] <- x$alpha.GammaX.iIB %*% x$dLambda.dtheta[[iName]]
        }else if(iType == "B"){
            x$dmu.dtheta[[iName]] <- x$alpha.GammaX.iIB %*% x$dB.dtheta[[iName]] %*% x$iIB.Lambda
        }
        
        ## *** variance-covariance matrix
        if(iType %in% "Psi_var"){
            x$dOmega.dtheta[[iName]] <-  t(x$iIB.Lambda) %*% x$dPsi.dtheta[[iName]] %*% x$iIB.Lambda
        }else if(iType %in% "Psi_cov"){
            x$dOmega.dtheta[[iName]] <-  t(x$iIB.Lambda) %*% x$dPsi.dtheta[[iName]] %*% x$iIB.Lambda
        }else if(iType == "Lambda"){
            x$dOmega.dtheta[[iName]] <- x$tLambda.tiIB.Psi.iIB %*% x$dLambda.dtheta[[iName]]
            x$dOmega.dtheta[[iName]] <- x$dOmega.dtheta[[iName]] + t(x$dOmega.dtheta[[iName]])
        }else if(iType == "B"){
            x$dOmega.dtheta[[iName]] <- x$tLambda.tiIB.Psi.iIB %*% x$dB.dtheta[[iName]] %*% x$iIB.Lambda
            x$dOmega.dtheta[[iName]] <- x$dOmega.dtheta[[iName]] + t(x$dOmega.dtheta[[iName]])
        }

### ** debug
        ## Omega <- moments(x, p = pars(x), conditional=TRUE, data = data)$C
        ## iOmega <- solve(Omega)
        ## tiOmega.epsilon <- t(iOmega %*% t(residuals(x)))
        
        ## if(!is.null(dmu.dtheta[[iP]])){
        ##     iScore <- (dmu.dtheta[[iP]] * tiOmega.epsilon)
        ## }
        ## if(!is.null(dOmega.dtheta[[iP]])){
        ##     firstTerm <- - tr(iOmega %*% dOmega.dtheta[[iP]])            
        ##     secondTerm <- (tiOmega.epsilon %*% dOmega.dtheta[[iP]]) * tiOmega.epsilon
        ##     iScore <- as.double(firstTerm) + secondTerm            
        ## }
    }

    
### ** export
    return(x)
}

## * .calcLeverage
.calcLeverage <- function(n.group, table.group, ls.indexGroup,
                          dmu.dtheta, vcov.param,
                          Omega, iOmega, Omega_chol,
                          power, as.clubSandwich){

    ls.leverage <- lapply(1:n.group, function(iG){ # iG <- 1

### ** Prepare
        n.tempo <- table.group[[iG]]
        Id.tempo <- diag(1,nrow=n.tempo,ncol=n.tempo)

        if(!is.null(ls.indexGroup)){
            index.tempo <- ls.indexGroup[[iG]]            
            iOmega.tempo <- iOmega[index.tempo,index.tempo,drop=FALSE]
            dmu.dtheta.tempo <- dmu.dtheta[index.tempo,,drop=FALSE]
            if(power != 1 && as.clubSandwich){
                Omega_chol.tempo <- Omega_chol[index.tempo,index.tempo,drop=FALSE]
                Omega.tempo <- Omega[index.tempo,index.tempo,drop=FALSE] 
            }
        }else{
            iOmega.tempo <- iOmega            
            dmu.dtheta.tempo <- dmu.dtheta[(n.tempo*(iG-1)+1):(n.tempo*iG),,drop=FALSE]
            if(power != 1 && as.clubSandwich){
                Omega_chol.tempo <- Omega_chol
                Omega.tempo <- Omega
            }
        }
### ** Compute IH
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




#----------------------------------------------------------------------
### score2.R ends here

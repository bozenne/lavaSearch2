### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: okt 25 2017 (12:06) 
##           By: Brice Ozenne
##     Update #: 1137
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
score2.gls <- function(x, cluster, p = NULL, data = NULL, indiv = TRUE, adjust.residuals = TRUE, power = 1/2,
                       as.clubSandwich = TRUE, ...){

    ## ** Normalize arguments
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
    
    x.coef <- coef(x)
    if(!is.null(p) && any(abs(p-x.coef)>1e-10)){
        stop("can only compute the score at the estimated model parameters \n",
             "consider setting argument \'p\' to NULL \n")
    }

    ## ** identify groups
    table.group <- table(vec.group)
    Ugroup <- unique(vec.group)
    n.group <- length(Ugroup)

    ls.indexGroup <- lapply(1:n.group, function(iG){
        which(vec.group==iG)
    })

    
    ## ** Reconstruct covariance matrix
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
    ## check that Sigma_beta = (t(X) Sigma_epsilon^{-1} X)^{-1}
    ## vcov(x) / solve(t(X) %*% solve(Matrix::bdiag(ls.Omega)) %*% X)

    ## ** Reconstruct variance covariance matrix (parameters)
    factor <- (x$dims$N - (x$method == "REML") * x$dims$p)/(x$dims$N-x$dims$p)
    vcov.param <- vcov(x) / factor

    ## ** Compute score
    score0 <- .score2_nlme(name.param = names(x.coef),
                           dmu.dtheta = X, epsilon = residuals(x, type = "response", level = 0),
                           vcov.param = vcov.param, ls.Omega = ls.Omega,
                           ls.indexGroup = ls.indexGroup, n.group = n.group, table.group = table.group,
                           indiv = indiv, adjust.residuals = adjust.residuals, power = power, as.clubSandwich = as.clubSandwich)
    
    ## ** Export
    return(score0)
}


## * score2.lme
score2.lme <- function(x, p = NULL, data = NULL, indiv = TRUE, adjust.residuals = TRUE, power = 1/2,
                       as.clubSandwich = TRUE, ...){

    ## ** Normalize arguments
    if(!identical(class(x),"lme")){
        wrongClass <- paste(setdiff(class(x),"lme"), collapse = " ")
        stop("iid2 is not available for ",wrongClass," objects \n")
    }
    if(!is.null(x$modelStruct$corStruct)){
        stop("cannot handle lme objects when corStruct is not null \n")
    }
    if(NCOL(x$groups)>1){
        stop("cannot handle lme objects with more than one random effect \n")
    }
    vec.group <- as.numeric(x$groups[,1])
    if(all(diff(vec.group)>=0) == FALSE && all(diff(vec.group<=0)) == FALSE ){
        stop("the grouping variable for the random effect must be sorted before fitting the model \n")
    }
    x.coef <- fixef(x)
    if(!is.null(p) && any(abs(p-x.coef)>1e-10)){
        stop("can only compute the score at the estimated model parameters \n",
             "consider setting argument \'p\' to NULL \n")
    }
    
    if(is.null(data)){
        data <- getData(x)        
    }
    X <- model.matrix(formula(x), data)
    
    ## ** Prepare
    table.group <- table(vec.group)
    Ugroup <- unique(vec.group)
    n.group <- length(Ugroup)

    ls.indexGroup <- lapply(1:n.group, function(iG){
        which(vec.group==iG)
    })

    ## ** Reconstruct covariance matrix
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

    ## check that Sigma_beta = (t(X) Sigma_epsilon^{-1} X)^{-1}
    ## vcov(x) - solve(t(X) %*% solve(Matrix::bdiag(ls.Omega)) %*% X)

    ## ** Compute score
    score0 <- .score2_nlme(name.param = names(x.coef),
                           dmu.dtheta = X, epsilon = residuals(x, type = "response", level = 0),
                           vcov.param = vcov(x), ls.Omega = ls.Omega,
                           ls.indexGroup = ls.indexGroup, n.group = n.group, table.group = table.group,
                           indiv = indiv, adjust.residuals = adjust.residuals, power = power, as.clubSandwich = as.clubSandwich)
    
    ## ** Export
    return(score0)
}

## * score2.lvmfit
score2.lvmfit <- function(x, p = NULL, data = NULL, indiv = TRUE,
                          adjust.residuals = TRUE, power = 1/2,
                          as.clubSandwich = TRUE, ...){

### ** normalize arguments
    if(!identical(class(x),"lvmfit")){
        wrongClass <- paste(setdiff(class(x),"lvmfit"), collapse = " ")
        stop("score2 is not available for ",wrongClass," objects \n")
    }
    if(!is.null(x$model0$attributes$type)){
        stop("score2 is only available for latent variable models involving gaussian variables \n")
    }
    
    if(is.null(p)){
        p <- pars(x)
    }else{
        if(is.null(names(p))){
            stop("\'p\' must be names \n")
        }
        ref.tempo <- attr(coefType(x),"reference")
        expected.names <- names(ref.tempo)[ref.tempo==FALSE] 
        if(any(expected.names %in% names(p) == FALSE)){ ## used in  .calcDtheta
            stop("argument \'p\' does not contain all the necessary parameters \n")
        }
    }

    if(is.null(data)){
        data <- model.frame(x)
    }
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }
    
### ** prepare elements

    ##
    name.param <- names(p)
    n.param <- length(p)

    n <- NROW(data)

    name.endogenous <- endogenous(x)
    n.endogenous <- length(name.endogenous)

    ## Omega
    Omega <- moments(x, p = p, conditional=TRUE, data = data)$C
    Omega_chol <- chol(Omega)
    iOmega <- chol2inv(Omega_chol) ## faster compared to solve
    
### ** compute derivatives
    ls.calc <- .calcDtheta(x, data = data,
                           param = p, n.param = n.param,
                           dmu = TRUE, dOmega = TRUE)
                                        # names(ls.calc$dmu.dtheta)
### ** residuals
    epsilon <- as.matrix(residuals(x, p = p))
    
    if(adjust.residuals){

        vcov.param <- vcov(x)
        ## *** gather derivatives
        dmu.dtheta <- sapply(colnames(vcov.param), function(iP){
            if(is.null(ls.calc$dmu.dtheta[[iP]])){
                return(rep(0, times = n * n.endogenous))
            }else{
                return(as.vector(t(ls.calc$dmu.dtheta[[iP]])))
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
   
### ** compute score
    n.score <- indiv*n+(1-indiv)*1 ## n if indiv==TRUE and 1 if indiv==FALSE
    out.score <- matrix(NA, nrow = n.score, ncol = n.param, dimnames = list(NULL,name.param))
    fct.sum <- switch(as.character(as.double(indiv)),
                      "1" = "rowSums",
                      "0" = "sum")

    ## tiOmega.epsilon <- t(iOmega %*% t(epsilon))
    tiOmega.epsilon <- epsilon %*% iOmega

    for(iP in 1:n.param){ # iP <- 1
        iScore <- rep(0,n.score)
        if(!is.null(ls.calc$dmu.dtheta[[iP]])){
            iScore <- do.call(fct.sum, args = list(ls.calc$dmu.dtheta[[iP]] * tiOmega.epsilon))
        }

        if(!is.null(ls.calc$dOmega.dtheta[[iP]])){
            firstTerm <- - (indiv*1+(1-indiv)*n)/2 * tr(iOmega %*% ls.calc$dOmega.dtheta[[iP]])            
            secondTerm <- 1/2 * do.call(fct.sum, args = list((tiOmega.epsilon %*% ls.calc$dOmega.dtheta[[iP]]) * tiOmega.epsilon))
            iScore <- iScore + as.double(firstTerm) + secondTerm            
        }
        out.score[,iP] <- iScore       
    }

### ** export
    return(out.score)

}

## * .score2_nlme
.score2_nlme <- function(name.param,
                         dmu.dtheta, epsilon, vcov.param, ls.Omega,  
                         indiv, adjust.residuals, power, as.clubSandwich,
                         ls.indexGroup, n.group, table.group){

    Omega <- Matrix::bdiag(ls.Omega)
    iOmega <- Matrix::bdiag(lapply(ls.Omega,solve))
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
        score.tempo <- t(dmu.dtheta[index.tempo,,drop=FALSE]) %*% iOmega[index.tempo,index.tempo,drop=FALSE] %*% epsilon[index.tempo]
        return(as.double(score.tempo))
    }))
   
    ## ** export
    if(indiv==FALSE){
        score0 <- rbind(colSums(score0))
    }
    colnames(score0) <- name.param

    return(score0)
   
}

#----------------------------------------------------------------------
### score2.R ends here

## * .calcDetheta
.calcDtheta <- function(x, data, param, n.param, 
                        dmu, dOmega){

    n <- NROW(data)

    allType <- coefType(x, detailed = TRUE)
    type <- allType[attr(allType,"reference")==FALSE]

    name.param <- names(param)

    name.endogenous <- endogenous(x)
    n.endogenous <- length(name.endogenous)

    name.latent <- latent(x)
    n.latent <- length(name.latent)
    
    name.allType <- names(allType)
    link.allType <- initVarLinks(name.allType)

    name.type <- names(type)
    link.type <- initVarLinks(name.type)
    names(link.type$var1) <- name.type
    names(link.type$var2) <- name.type
    
### ** precompute quantities
    
    ## *** coef
    allCoef <- coef(x,level=9)[,"Estimate"]
    allCoef[name.param] <- param  ## useful when the parameter values are defined by the user
    
    if(n.latent>0){
        ## *** Lambda
        Lambda <- matrix(0,nrow = n.latent, ncol = n.endogenous,
                           dimnames = list(name.latent,name.endogenous))
        index.Lambda <- which(allType=="Lambda")
        for(iLambda in index.Lambda){ # iLambda <- index.Lambda[1]
            Lambda[link.allType$var2[iLambda],link.allType$var1[iLambda]] <- allCoef[name.allType[iLambda]]
        }
        
        ## *** B
        B <- matrix(0,nrow = n.latent, ncol = n.latent,
                      dimnames = list(name.latent,name.latent))
        index.B <- which(allType=="B")
        for(iB in index.B){ # iB <- index.B[1]
            B[link.allType$var1[iB],link.allType$var2[iB]] <- allCoef[name.allType[iB]]
        }
        iB <- solve(diag(1,n.latent)-B)
        tiB.Lambda <- t(iB) %*% Lambda

        ## *** Psi
        Psi <- matrix(0,nrow = n.latent, ncol = n.latent,
                        dimnames = list(name.latent,name.latent))
        index.Psi <- which(allType == "Psi_var")
        if(length(index.Psi)>0){
            for(iPsi in index.Psi){ # iPsi <- index.Psi[1]
                Psi[link.allType$var1[iPsi],link.allType$var2[iPsi]] <- allCoef[name.allType[iPsi]]
            }
        }
        
        index.Psi <- which(allType == "Psi_cov")
        if(length(index.Psi)>0){
            for(iPsi in index.Psi){ # iPsi <- index.Psi[1]
                Psi[link.allType$var1[iPsi],link.allType$var2[iPsi]] <- allCoef[name.allType[iPsi]]
                Psi[link.allType$var2[iPsi],link.allType$var1[iPsi]] <- allCoef[name.allType[iPsi]]
            }
        }
        Psi.tiB <- Psi %*% t(iB)
        tLambda.iB.Psi.tiB <- t(tiB.Lambda) %*% Psi.tiB

        ## *** alpha Gamma X
        alpha.GammaX <- matrix(0,nrow = n, ncol = n.latent,
                                 dimnames = list(NULL,name.latent))

        for(iLatent in 1:n.latent){ # iLatent <- 1
            index.alpha <- intersect(which(link.allType$var1==name.latent[iLatent]),
                                     which(allType=="alpha"))
            if(length(index.alpha)>0){
                alpha.GammaX[,iLatent] <- allCoef[index.alpha]
            }
            
            index.Gamma <- intersect(which(link.allType$var1==name.latent[iLatent]),
                                     which(allType=="Gamma"))
            if(length(index.Gamma)>0){
                alpha.GammaX[,iLatent] <- alpha.GammaX[,iLatent] + data[,link.allType$var2[index.Gamma],drop=FALSE] %*% allCoef[index.Gamma]
            }
           
        }
        alpha.GammaX.iB <- alpha.GammaX %*% iB
    }
    
### ** partial derivative
    dmu.dtheta <- vector(mode = "list", length = n.param)
    dOmega.dtheta <- vector(mode = "list", length = n.param)
    names(dmu.dtheta) <- name.param
    names(dOmega.dtheta) <- name.param
    
    for(iP in 1:n.param){ # iP <- 1
        iType <- type[name.param[iP]]
        iVar1 <- link.type$var1[name.param[iP]]
        iVar2 <- link.type$var2[name.param[iP]]

        ## *** linear predictor
        if(dmu){
            if(iType == "nu"){                
                dmu.dtheta[[iP]] <- matrix(0,nrow=n,ncol=n.endogenous)
                dmu.dtheta[[iP]][,name.endogenous == iVar1] <- rep(1,n)
            }else if(iType == "K"){
                dmu.dtheta[[iP]] <- matrix(0,nrow=n,ncol=n.endogenous)
                dmu.dtheta[[iP]][,name.endogenous == iVar1] <- data[,iVar2]
            }else if(iType == "alpha"){
                dmu.dtheta[[iP]] <- matrix(as.numeric(tiB.Lambda[name.latent == iVar1,]),
                                           nrow = n, ncol = n.endogenous, byrow = TRUE)
            }else if(iType == "Gamma"){
                dmu.dtheta[[iP]] <- data[,iVar2] %o% tiB.Lambda[name.latent == iVar1,]
            }else if(iType == "Lambda"){
                J.tempo <- matrix(0,nrow = n.endogenous, ncol = n.latent)
                J.tempo[name.endogenous == iVar1,name.latent == iVar2] <- 1
                dmu.dtheta[[iP]] <- alpha.GammaX %*% t(iB) %*% t(J.tempo)
                ## alpha.GammaX.iB %*% t(J.tempo)
            }else if(iType == "B"){
                J.tempo <- matrix(0,nrow = n.latent, ncol = n.latent)
                J.tempo[name.latent == iVar1, name.latent == iVar2] <- 1
                dmu.dtheta[[iP]] <- alpha.GammaX %*% t(iB) %*% t(J.tempo) %*% tiB.Lambda
                ##dmu.dtheta[[iP]] <- - alpha.GammaX %*% iB %*% J.tempo %*% iB %*% t(Lambda)
            }
   
        }
        
        ## *** variance-covariance matrix
        if(dOmega){
            if(iType %in% c("Sigma_var")){
                dOmega.dtheta[[iP]] <- matrix(0,nrow=n.endogenous,ncol=n.endogenous)
                dOmega.dtheta[[iP]][name.endogenous == iVar1,name.endogenous == iVar2] <- 1
            }else if(iType == "Sigma_cov"){
                dOmega.dtheta[[iP]] <- matrix(0,nrow=n.endogenous,ncol=n.endogenous)
                dOmega.dtheta[[iP]][name.endogenous == iVar1,name.endogenous == iVar2] <- 1
                dOmega.dtheta[[iP]][name.endogenous == iVar2,name.endogenous == iVar1] <- 1
            }else if(iType %in% "Psi_var"){
                J.tempo <- matrix(0,nrow = n.latent, ncol = n.latent)
                J.tempo[name.latent == iVar1,name.latent == iVar2] <- 1
                dOmega.dtheta[[iP]] <-  t(tiB.Lambda) %*% J.tempo %*% tiB.Lambda
            }else if(iType %in% "Psi_cov"){
                J.tempo <- matrix(0,nrow = n.latent, ncol = n.latent)
                J.tempo[name.latent == iVar1,name.latent == iVar2] <- 1
                J.tempo[name.latent == iVar2,name.latent == iVar1] <- 1
                dOmega.dtheta[[iP]] <-  t(tiB.Lambda) %*% J.tempo %*% tiB.Lambda
            }else if(iType == "Lambda"){
                J.tempo <- matrix(0,nrow = n.endogenous, ncol = n.latent)
                J.tempo[name.endogenous == iVar1,name.latent == iVar2] <- 1
                dOmega.dtheta[[iP]] <- tLambda.iB.Psi.tiB %*% t(J.tempo)
                dOmega.dtheta[[iP]] <- dOmega.dtheta[[iP]] + t(dOmega.dtheta[[iP]])
            }else if(iType == "B"){
                J.tempo <- matrix(0,nrow = n.latent, ncol = n.latent)
                J.tempo[name.latent == iVar1, name.latent == iVar2] <- 1
                ## Omega - Lambda %*% iB %*% Psi %*% t(iB) %*% t(Lambda)
                ## dOmega.dtheta <-  Lambda %*% iB %*% Psi %*% t(iB) %*% t(J.tempo) %*% t(iB) %*% t(Lambda)
                dOmega.dtheta[[iP]] <- tLambda.iB.Psi.tiB %*% t(J.tempo) %*% tiB.Lambda
                dOmega.dtheta[[iP]] <- dOmega.dtheta[[iP]] + t(dOmega.dtheta[[iP]])
            }

    
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
    return(list(dmu.dtheta = dmu.dtheta,
                dOmega.dtheta = dOmega.dtheta))
}

## * .calcAdjustedResiduals
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




### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: okt 19 2017 (18:35) 
##           By: Brice Ozenne
##     Update #: 727
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
#' @param p [optional] vector of parameters.
#' @param data [optional] data set.
#' @param indiv Should the score relative to each observation be exported? Otherwise the total score (i.e. sum over all observations) will be exported.
#' @param adjust.residuals Should the leverage-adjusted residuals be used to compute the score? Otherwise the raw residuals will be used.
#' @param alpha the exponent used for adjusting the residuals: \eqn{e_adj = \frac{e}{(1-h_{ii})^\alpha}}.
#' \eqn{alpha=0.5} corresponds to leverage-adjustment (Kauermann and Carroll, 2001).
#' \eqn{alpha=1} approximates the leave-one-out jackknife variance estimator (Bell and McCaffrey, 2002) 
#' @param return.df Should the degree of freedom be exported?
#' 
#' @details A lvm can be written as a measurement model:
#' \deqn{Y_i = \nu + \eta_i \Lambda + X_i K + \epsilon_i}
#' and a structural model:
#' \deqn{\eta_i = \alpha + \eta_i B + X_i \Gamma + \zeta_i}
#' where \eqn{\Psi}   is the variance covariance matrix of the residuals \eqn{\zeta} \cr
#' and   \eqn{\Sigma} is the variance covariance matrix of the residuals \eqn{\epsilon}. \cr \cr
#'
#' \cr
#' 
#' The corresponding conditional mean is:
#' \deqn{
#' \mu_i(\theta) = E[Y_i|X_i] = \nu + (\alpha + X_i \Gamma) (1-B)^{-1} \Lambda + X_i K
#' }
#' \deqn{
#' \Omega(\theta) = Var[Y_i|X_i] = \Lambda^{t} (1-B)^{-t} \Psi (1-B)^{-1} \Lambda + \Sigma
#' }
#'
#' Therefore \eqn{\nu}, \eqn{K}, \eqn{\alpha}, \eqn{\Gamma} are pure mean parameters, \eqn{\Psi}, \eqn{\Sigma} pure variance parameters, \cr
#' and \eqn{\Lambda}, \eqn{B} are both mean and variance parameters. \cr
#'
#' \cr
#' 
#' The log-likelihood can written:
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
#' @references
#' Bell, R. M., & McCaffrey, D. F. Bias reduction in standard errors for linear regression with multi-stage samples. Survey Methodology, 28(2), 169-181 (2002). \cr
#' Kauermann G. and Carroll R. J. A note on the efficiency of sandwich covariance matrix estimation. Journal of the American Statistical Association. Vol. 96, No. 456 (2001).
#' @export
`score2` <-
  function(x, ...) UseMethod("score2")

## * score2.lvmfit
score2.lvmfit <- function(x, p = NULL, data = NULL, power = 1/2,
                          indiv = TRUE, adjust.residuals = TRUE, return.df = TRUE, ...){

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
    iOmega <- solve(Omega)
    
### ** compute derivatives
    ls.calc <- .calcDtheta(x, data = data,
                           param = p, n.param = n.param,
                           dmu = TRUE, dOmega = TRUE)
    # names(ls.calc$dmu.dtheta)

### ** residuals
    epsilon <- residuals(x, p = p)
    
    if(adjust.residuals){

        vcov.param <- vcov(x)
        ## *** gather derivatives
        ls.M <- lapply(ls.calc$dmu.dtheta[colnames(vcov.param)], function(iD){
            if(is.null(iD)){
                return(matrix(0,nrow = n,ncol=n.endogenous))
            }else{
                return(iD)
            }
        })

        ## *** compute leverage values
        if(n.endogenous==1){
            dmu.dtheta2 <- do.call(cbind, args = ls.M)
            leverage.adj <- 1 - diag(iOmega) * rowSums( (dmu.dtheta2 %*% vcov.param) * dmu.dtheta2)
        }else{
            dmu.dtheta2 <- do.call(abind::abind, args = list(ls.M, along = 3))

            ls.diagH <- lapply(1:n, function(iI){ # iI <- 1
                ## rowSums((dmu.dtheta2[iI,,] %*%  vcov.param) * dmu.dtheta2[iI,,]) * diag(iOmega)
                dmu.dtheta2[iI,,] %*%  vcov.param %*% t(dmu.dtheta2[iI,,]) %*% diag(iOmega)
            })
            
            leverage.adj <- 1 - t(do.call(cbind,ls.diagH))

        }
        ## *** normalize residuals
        epsilon <- epsilon/leverage.adj^(power)
        
    }
    #tiOmega.epsilon <- t(iOmega %*% t(epsilon))
    tiOmega.epsilon <- epsilon %*% iOmega
   
### ** compute score
    n.score <- indiv*n+(1-indiv)*1 ## n if indiv==TRUE and 1 if indiv==FALSE
    out.score <- matrix(NA, nrow = n.score, ncol = n.param, dimnames = list(NULL,name.param))
    fct.sum <- switch(as.character(as.double(indiv)),
                      "1" = "rowSums",
                      "0" = "sum")
    
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

### ** degree of freedom
    if(return.df){
        df.adj <- df.residual(x, adjust.residuals = adjust.residuals, 
                              dmu.dtheta2 = dmu.dtheta2, leverage.adj = leverage.adj, Omega = Omega, iOmega = iOmega,
                              p = p)
        attr(out.score,"df") <- df.adj
    }

### ** export
    ## out.score - score(x, indiv = TRUE)[,colnames(out.score)]
    
    return(out.score)

}

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


#----------------------------------------------------------------------
### score2.R ends here

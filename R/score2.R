### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: okt 13 2017 (16:10) 
##           By: Brice Ozenne
##     Update #: 257
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Compute the score directly from the gaussian density
#' @description Compute the score directly from the gaussian density
#'
#' 
#' @details A lvm can be written as a measurement model:
#' \deqn{Y_i = \nu + \Lambda \eta_i + K X_i + \epsilon_i}
#' and a structural model:
#' \deqn{\eta_i = \alpha + B \eta_i + \Gamma X_i + \zeta_i}
#' where \eqn{\Psi}   is the variance covariance matrix of the residuals \eqn{\zeta} \cr
#' and   \eqn{\Sigma} is the variance covariance matrix of the residuals \eqn{\epsilon}. \cr \cr
#'
#' \cr
#' 
#' The corresponding conditional mean is:
#' \deqn{
#' \mu_i(\theta) = E[Y_i|X_i] = \nu + \Lambda (1-B)^{-1} \alpha + (\Lambda (1-B)^{-1} \Gamma + K) X_i 
#' }
#' \deqn{
#' \Omega(\theta) = Var[Y_i|X_i] = \Lambda (1-B)^{-1} \Psi (1-B)^{-t} \Lambda^t + \Sigma
#' }
#'
#' Therefore \eqn{\nu}, \eqn{K}, \eqn{\alpha}, \eqn{\Gamma} are pure mean parameters, \eqn{\Psi}, \eqn{\Sigma} pure variance parameters, \cr
#' and \eqn{\Lambda}, \eqn{B} pure variance parameters are both mean and variance parameters. \cr
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
#' \deqn{\frac{\partial \mu_i(\theta)}{\partial \nu} = 1}
#' \deqn{\frac{\partial \mu_i(\theta)}{\partial \K} = X_i}
#' \deqn{\frac{\partial \mu_i(\theta)}{\partial \alpha} = \Lambda(1-B)^{-1}}
#' \deqn{\frac{\partial \mu_i(\theta)}{\partial \Gamma} = \Lambda(1-B)^{-1}X_i}
#' \deqn{\frac{\partial \mu_i(\theta)}{\partial \Lambda} = (1-B)^{-1}(\alpha + \Gamma X_i)}
#' \deqn{\frac{\partial \mu_i(\theta)}{\partial B} = -\Lambda(1-B)^{-1}(1-B)^{-1}(\alpha + \Gamma X_i)}
#' and:
#' \deqn{\frac{\partial \Omega(\theta)}{\partial \psi} = \Lambda (1-B)^{-1} \delta_{\psi \in \Psi} (1-B)^{-t} \Lambda^t}
#' \deqn{\frac{\partial \Omega(\theta)}{\partial \sigma} = \delta_{\sigma \in \Sigma}}
#' \deqn{\frac{\partial \Omega(\theta)}{\partial \lambda} = \delta_{\lambda \in \Lambda} (1-B)^{-1} \Psi (1-B)^{-t} \Lambda^t + \Lambda (1-B)^{-1} \Psi (1-B)^{-t} \delta_{\lambda \in \Lambda}}
#' \deqn{\frac{\partial \Omega(\theta)}{\partial B} = - \Lambda (1-B)^{-1} \delta_{b \in B} (1-B)^{-1} \Psi (1-B)^{-t} \Lambda^t - \Lambda (1-B)^{-1} \Psi (1-B)^{-1,t} \delta_{b \in B} (1-B)^{-1,t}}
#' @examples
#'
#' @export
`score2` <-
  function(x, ...) UseMethod("score2")

score2.lvmfit <- function(x, mu = NULL, Omega = NULL, data = NULL, adjust.residuals = TRUE, ...){

### ** normalize arguments
    if(!identical(class(x),"lvmfit")){
        wrongClass <- paste(setdiff(class(x),"lvmfit"), collapse = " ")
        stop("score2 is not available for ",wrongClass," objects \n")
    }
    allType <- coefType(x, detailed = TRUE)
    type <- allType[attr(allType,"reference")==FALSE]

    if(is.null(data)){
        data <- model.frame(x)
    }
    data <- as.matrix(data)
    if(is.null(mu) || is.null(Omega)){
        mom.th <- moments(x, p = pars(x), conditional=TRUE, data = data)
        if(is.null(mu)){
            mu <- mom.th$xi
        }
        if(is.null(Omega)){
            Omega <- mom.th$C
        }    
    }

### ** tag elements
    p <- length(type)
    n <- NROW(data)
    name.endogenous <- endogenous(x)
    n.endogenous <- length(name.endogenous)
    name.latent <- latent(x)
    n.latent <- length(name.latent)
    
    name.allType <- names(allType)
    link.allType <- initVarLinks(name.allType)

    name.type <- names(type)
    link.type <- initVarLinks(name.type)

### ** precompute quantities

    ## *** coef
    allCoef <- coef(x,level=9)[,"Estimate"]
    
    ## *** residuals
    epsilon <- residuals(x)

    if(adjust.residuals){
        
    }
    ## *** iOmega
    iOmega <- solve(Omega)
    epsilon.iOmega <- epsilon %*% iOmega
    
    if(n.latent>0){
        ## *** Lambda
        M.Lambda <- matrix(0,nrow = n.endogenous, ncol = n.latent,
                           dimnames = list(name.endogenous,name.latent))
        index.Lambda <- which(allType=="Lambda")
        for(iLambda in index.Lambda){ # iLambda <- index.Lambda[1]
            M.Lambda[link.allType$var1[iLambda],link.allType$var2[iLambda]] <- allCoef[name.allType[iLambda]]
        }
        
        ## *** B
        M.B <- matrix(0,nrow = n.latent, ncol = n.latent,
                      dimnames = list(name.latent,name.latent))
        index.B <- which(allType=="B")
        for(iB in index.B){ # iB <- index.B[1]
            M.B[link.allType$var1[iB],link.allType$var2[iB]] <- allCoef[name.allType[iB]]
        }
        iB <- solve(diag(1,n.latent)-M.B)
        Lambda.iB <- M.Lambda %*% iB
        
        ## *** Psi
        M.Psi <- matrix(0,nrow = n.latent, ncol = n.latent,
                        dimnames = list(name.latent,name.latent))
        index.Psi <- which(allType %in% c("Psi_var","Psi_cov"))
        for(iPsi in index.Psi){ # iPsi <- index.Psi[1]
            M.Psi[link.allType$var1[iPsi],link.allType$var2[iPsi]] <- allCoef[name.allType[iPsi]]
        }
        Psi.iB <- M.Psi %*% iB
        iB.Lambda.iB.Psi <- Lambda.iB %*% Psi.iB

        ## *** alpha Gamma X
        M.alpha.GammaX <- matrix(0,nrow = n, ncol = n.latent,
                                 dimnames = list(NULL,name.latent))

        for(iLatent in 1:n.latent){ # iLatent <- 1
            index.alpha <- intersect(which(link.allType$var1==name.latent[iLatent]),
                                     which(allType=="alpha"))
            if(length(index.alpha)>0){
                M.alpha.GammaX[,iLatent] <- allCoef[index.alpha]
            }
            
            index.Gamma <- intersect(which(link.allType$var1==name.latent[iLatent]),
                                     which(allType=="Gamma"))
            if(length(index.Gamma)>0){
                M.alpha.GammaX[,iLatent] <- M.alpha.GammaX[,iLatent] + data[,link.allType$var2[index.Gamma]] %*% allCoef[index.Gamma]
            }
           
        }

        iB.alpha.GammaX <- M.alpha.GammaX %*% iB
    }
    
    

### ** compute derivatives
    out.score <- matrix(NA,nrow = n, ncol = p, dimnames = list(NULL,name.type))

    for(iP in 1:p){ # iP <- 1
        iType <- type[iP]
        iName <- name.type[iP]
        iVar <- lapply(link.type,"[",iP)
        print(iName)
        
        ## *** partial derivative of the linear predictor
        if(iType == "nu"){
            dmu.dtheta <- matrix(0,nrow=n,ncol=n.endogenous)
            dmu.dtheta[,name.endogenous == iVar$var1] <- rep(1,n)
            add.mean <- TRUE
        }else if(iType == "alpha"){
            dmu.dtheta <- matrix(as.numeric(Lambda.iB[,name.latent == iVar$var1]),
                                 nrow = n, ncol = n.endogenous, byrow = TRUE)
            add.mean <- TRUE
        }else if(iType == "Gamma"){
            dmu.dtheta <- data[,iVar$var2] %o% Lambda.iB[,iVar$var1]
            add.mean <- TRUE
        }else if(iType == "K"){
            dmu.dtheta <- matrix(0,nrow=n,ncol=n.endogenous)
            dmu.dtheta[,name.endogenous == iVar$var1] <- data[,iVar$var2]
            add.mean <- TRUE
        }else if(iType == "Lambda"){
            browser()
            dmu.dtheta <- matrix(0,nrow=n,ncol=n.endogenous)
            dmu.dtheta[,name.endogenous==iVar$var1] <- iB.alpha.GammaX %*% (name.latent == iVar$var2)
            add.mean <- TRUE
        }else if(iType == "B"){
            browser()            
            ## dmu.dtheta <- matrix(0,nrow=n,ncol=n.endogenous)
            ## dmu.dtheta[,name.endogenous==iVar$var1] <- Lambda.iB[,name.latent == iVar$var2] iB.alpha.GammaX %*% (name.latent == iVar$var2)
            ## dmu.dtheta <- matrix(0,nrow=n,ncol=n.endogenous)
            add.mean <- TRUE
        }else{
            add.mean <- FALSE
        }

        ## *** partial derivative of the covariance
        dOmega.dtheta <- matrix(0,nrow=n.endogenous,ncol=n.endogenous)
        if(iType %in% c("Sigma_var")){
            dOmega.dtheta[name.endogenous == iVar$var1,name.endogenous == iVar$var2] <- 1
            add.cov <- TRUE
        }else if(iType == "Sigma_cov"){
            dOmega.dtheta[name.endogenous == iVar$var1,name.endogenous == iVar$var2] <- 1
            dOmega.dtheta[name.endogenous == iVar$var2,name.endogenous == iVar$var1] <- 1
            add.cov <- TRUE
        }else if(iType %in% "Psi_var"){
            J.tempo <- matrix(0,nrow = n.latent, ncol = n.latent)
            J.tempo[name.latent == iVar$var1,name.latent == iVar$var2] <- 1
            dOmega.dtheta <-  Lambda.iB %*% J.tempo %*% t(Lambda.iB)
            add.cov <- TRUE
        }else if(iType %in% "Psi_cov"){
            J.tempo <- matrix(0,nrow = n.latent, ncol = n.latent)
            J.tempo[name.latent == iVar$var1,name.latent == iVar$var2] <- 1
            J.tempo[name.latent == iVar$var2,name.latent == iVar$var1] <- 1
            dOmega.dtheta <-  Lambda.iB %*% J.tempo %*% t(Lambda.iB)
            add.cov <- TRUE
        }else if(iType == "Lambda"){
            browser()
            dOmega.dtheta <- iB.Lambda.iB.Psi[,name.latent == iVar$var2,drop=FALSE] %*% (name.endogenous == iVar$var1)
            dOmega.dtheta <- dOmega.dtheta + t(dOmega.dtheta)
                                        #            dOmega.dtheta <- iB.Lambda.iB.Psi %*% (name.endogenous == iVar$var1) + t(iB.Lambda.iB.Psi %*% (name.endogenous == iVar$var1))
            add.cov <- TRUE
        }else if(iType == "B"){
            browser()
            add.cov <- TRUE
        }else{
            add.cov <- FALSE
        }
        
        
      
        ## *** add contributions
        iScore <- rep(0,n)

        if(add.mean){
            iScore <- rowSums(dmu.dtheta * epsilon.iOmega)
        }

        
        if(add.cov){
            firstTerm <- - 1/2 * tr(iOmega %*% dOmega.dtheta)
            for(iObs in 1:n){ # iObs <- 1
                iScore[iObs] <-  iScore[iObs] + firstTerm + 1/2 * t(epsilon.iOmega[iObs,]) %*% dOmega.dtheta %*% epsilon.iOmega[iObs,]
            }

        }
        
        out.score[,iName] <- iScore       
    }

    ## out.score - score(x, indiv = TRUE)[,colnames(out.score)]
    return(out.score)

}

#----------------------------------------------------------------------
### score2.R ends here

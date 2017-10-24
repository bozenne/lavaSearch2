### df.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 23 2017 (12:27) 
## Version: 
## last-updated: okt 24 2017 (10:21) 
##           By: Brice Ozenne
##     Update #: 106
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Approximate degree of freedoms of a model
#' @description Approximate the degree of freedoms of a model
#' @name df.residual
#' 
#' @param object model object
#' @param adjust.residuals should leverage-adjusted residuals be used? 
#' @param power the type of leverage-adjusted residuals to be used. See the details section.
#' @param dmu.dtheta2 
#' @param leverage.adj 
#' @param Omega 
#' @param iOmega 
#' @param data 
#' @param p 
#' @param ... additional arguments
#'
#' @details
#' Leverage-adjusted residuals have been shown to improve the coverage of robust standard errors in small samples.
#' 
#' Since the calculation of the degrees of freedom involves the use of the residuals,
#' the same small sample correction should be use to compute the degrees of freedom and the robust standard errors.
#'
#' 
#' 
#' 
#' @examples
#' n <- 20
#' 
#' set.seed(10)
#' m <- lvm(Y~X1+X2+X3+X4)
#' d <- lava::sim(m, n)
#' e <- estimate(m, data = d)
#' df.residual(e)
#'
#' 
#' @rdname df.residual
#' @method df.residual lvmfit
#' @export
df.residual.lvmfit <- function(object, power = 1/2, adjust.residuals=TRUE,
                               dmu.dtheta2=NULL, leverage.adj=NULL, Omega=NULL, iOmega=NULL,
                               data=NULL,p=NULL,
                               ...) {

### ** normalize arguments
    if(!identical(class(object),"lvmfit")){
        wrongClass <- paste(setdiff(class(object),"lvmfit"), collapse = " ")
        stop("score2 is not available for ",wrongClass," objects \n")
    }

    if(is.null(p)){
        p <- pars(object)
    }else{
        if(is.null(names(p))){
            stop("\'p\' must be names \n")
        }
        ref.tempo <- attr(coefType(object),"reference")
        expected.names <- names(ref.tempo)[ref.tempo==FALSE] 
        if(any(expected.names %in% names(p) == FALSE)){ ## used in  .calcDtheta
            stop("argument \'p\' does not contain all the necessary parameters \n")
        }
    }
### ** precompute quantities
    name.endogenous <- endogenous(object)
    n.endogenous <- length(name.endogenous)

    if(is.null(data)){
        data <- model.frame(object)
    }
    if(!is.matrix(data)){
        data <- as.matrix(data)
    }
    
    if(is.null(Omega)){
        Omega <- moments(object, p = p, conditional=TRUE, data = data)$C
    }

    if(is.null(iOmega)){
        iOmega <- solve(Omega)
    }

    vcov.param <- vcov(object)
    if(adjust.residuals && (is.null(dmu.dtheta2) || is.null(leverage.adj))){
        
        ## derivative of the linear predictor
        ls.calc <- .calcDtheta(object, data = data, param = p, n.param = length(p), 
                               dmu = TRUE, dOmega = FALSE)

        ## gather derivatives
        ls.M <- lapply(ls.calc$dmu.dtheta[colnames(vcov.param)], function(iD){
            if(is.null(iD)){
                return(matrix(0,nrow = n,ncol=n.endogenous))
            }else{
                return(iD)
            }
        })
        ## compute leverage values
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
    }
    

    

### ** compute degree of freedom
    if(adjust.residuals){
        if(n.endogenous==1){
            df.adj <- sum(leverage.adj)
            
        }else{
            eigen.Omega <- eigen(Omega)
            root.Omega <- eigen.Omega$vectors %*% diag(sqrt(eigen.Omega$values), nrow = n.endogenous, ncol = n.endogenous) %*% t(eigen.Omega$vectors)
            ## root.Omega %*% root.Omega - Omega

            dmu.dtheta <- apply(dmu.dtheta2, 3, function(iP){
              as.vector(t(iP))  
            })
            diag.iOmega <- Matrix::bdiag(lapply(1:n, function(i){iOmega}))
                
            H <- dmu.dtheta %*% vcov.param %*% t(dmu.dtheta) %*% diag.iOmega
            vec.eigen <- eigen(H)$values
            df.adj <- sum(unlist(vec.eigen))^2/sum(unlist(vec.eigen)^2)
            
            print(colSums(leverage.adj))
        }
    }else{
        df.adj <- n-1
    }


### ** export
    return(df.adj)
}

#' @rdname df.residual
#' @method df.residual coxph
#' @export
df.residual.coxph <- function(object, ...) {
    n <- riskRegression::coxN(object)    
    object.coef <- coef(object)
    p.effective <- length(object.coef)
    
    return(n-p.effective)
}
#----------------------------------------------------------------------
### df.R ends here

### df.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 23 2017 (12:27) 
## Version: 
## last-updated: okt 24 2017 (17:13) 
##           By: Brice Ozenne
##     Update #: 113
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
df.residual.lvmfit <- function(object, power=1/2, adjust.residuals=TRUE, as.clubSandwich=TRUE,
                               dmu.dtheta = NULL, vcov.param, ls.leverage.adj=NULL,
                               Omega=NULL, iOmega=NULL, Omega_chol=NULL,
                               ...) {

### ** Normalize arguments
    if(!identical(class(x),"lvmfit")){
        wrongClass <- paste(setdiff(class(x),"lvmfit"), collapse = " ")
        stop("score2 is not available for ",wrongClass," objects \n")
    }
    if(!is.null(x$model0$attributes$type)){
        stop("score2 is only available for latent variable models involving gaussian variables \n")
    }

    
### ** precompute quantities
    p <- pars(x)
    name.param <- names(p)
    n.param <- length(p)

    data <- as.matrix(model.frame(x))
    n <- NROW(data)

    name.endogenous <- endogenous(object)
    n.endogenous <- length(name.endogenous)

    vcov.param <- vcov(x)

    ## Omega
    Omega <- moments(x, p = p, conditional=TRUE, data = data)$C
    Omega_chol <- chol(Omega)
    iOmega <- chol2inv(Omega_chol) ## faster compared to solve
    
    ## *** compute derivatives
    ls.calc <- .calcDtheta(x, data = data,
                           param = p, n.param = n.param,
                           dmu = TRUE, dOmega = TRUE)

    ## *** gather derivatives
    dmu.dtheta <- sapply(colnames(vcov.param), function(iP){
        if(is.null(ls.calc$dmu.dtheta[[iP]])){
            return(rep(0, times = n * n.endogenous))
        }else{
            return(as.vector(t(ls.calc$dmu.dtheta[[iP]])))
        }
    })
    colnames(dmu.dtheta) <- colnames(vcov.param)

    ## *** compute leverage
    ls.leverage <- .calcLeverage(n.group = n, table.group = rep(n.endogenous, n), ls.indexGroup = NULL,
                                 dmu.dtheta = dmu.dtheta, vcov.param = vcov.param,
                                 Omega = Omega, iOmega = iOmega, Omega_chol = Omega_chol,
                                 power = power, as.clubSandwich = as.clubSandwich)


### ** compute df
    df.adj <- .calcDF(dmu.dtheta, vcov.param = vcov.param, ls.leverage = ls.leverage,
                      Omega = Omega, iOmega = iOmega, Omega_chol = Omega_chol)
    
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

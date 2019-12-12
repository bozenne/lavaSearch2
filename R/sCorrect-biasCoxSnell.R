### biasCoxSnell.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  2 2019 (10:20) 
## Version: 
## Last-Updated: dec 11 2019 (14:11) 
##           By: Brice Ozenne
##     Update #: 27
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

biasCoxSnell <- function(object){

    ## ** extract information from object
    param <- coef2(e)
    name.param <- names(param)
    n.param <- length(name.param)
    data <- model.frame(e)
    n <- nobs(e)

    ## *** identify missing values
    names(object)
    object$cluster
    
    ## *** pre-compute skeleton for conditional moments
    if(is.null(object$sCorrect$conditionalMoment)){
        object$sCorrect$conditionalMoment <- conditionalMoment(e, data = data,
                                                               first.order = TRUE, second.order = TRUE,
                                                               name.endogenous = endogenous(e),
                                                               name.latent = latent(e), usefit = FALSE)
    }

    
    ## ** update conditional moments
    moment.object <- conditionalMoment(e, data = data, param = param,
                                       first.order = TRUE, second.order = TRUE,
                                       name.endogenous = endogenous(e),
                                       name.latent = latent(e), usefit = TRUE)

    Omega <- moment.object$Omega
    dmu <- moment.object$dmu
    dOmega <- moment.object$dOmega
    d2mu <- moment.object$d2mu
    d2Omega <- moment.object$d2Omega

    ## ** update variance-covariance matrix of the parameters
    if(is.null(vcov)){
        grid.information <- .gridInformation(param = name.param,
                                             param.mean = attr(param, "mean.coef"),
                                             param.var = attr(param, "var.coef"))
        
        iVcov <- .vcov2(dmu = dmu,
                        dOmega = dOmega,
                        Omega = Omega,
                        n.corrected = n,
                        leverage = leverage, index.Omega = index.Omega, n.cluster = n.cluster,
                        grid.meanparam = grid.information$mean,
                        n.grid.meanparam = grid.information$n.mean,
                        grid.varparam = grid.information$var,
                        n.grid.varparam = grid.information$n.var,
                        name.param = name.param,
                        n.param = n.param,
                        attr.info = TRUE)
        iInfo <- .information2(dmu = dmu,
                               dOmega = dOmega,
                               Omega = Omega,
                               n.corrected = n.corrected,
                               leverage = leverage, index.Omega = index.Omega, n.cluster = n.cluster,
                               grid.meanparam = grid.meanparam,
                               n.grid.meanparam = n.grid.meanparam,
                               grid.varparam = grid.varparam,
                               n.grid.varparam = n.grid.varparam,
                               name.param = name.param,
                               n.param = n.param)
        iVcov.param <- try(chol2inv(chol(iInfo)), silent = TRUE)
        if(inherits(iVcov.param, "try-error")){
            iVcov.param <- solve(iInfo)
        }

    }
    

    ## 2*J(\theta3;\theta2,\theta1) + K(\theta3,\theta2,\theta1)    

    moment.object$value
}


######################################################################
### biasCoxSnell.R ends here

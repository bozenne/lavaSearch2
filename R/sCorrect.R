### sCorrect.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2018 (14:29) 
## Version: 
## Last-Updated: jan  7 2020 (13:58) 
##           By: Brice Ozenne
##     Update #: 1675
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - sCorrect
#' @title  Satterthwaite Correction and Small Sample Correction
#' @description Correct the bias of the ML estimate of the variance and compute the first derivative of the information matrix.
#' @name sCorrect
#'
#' @param object,x a \code{gls}, \code{lme}, or \code{lvm} object.
#' @param param [numeric vector, optional] the values of the parameters at which to perform the correction.
#' @param data [data.frame, optional] the dataset relative to which the correction should be performed.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' Only required for \code{gls} models with no correlation argument.
#' @param value [logical] value for the arguments \code{adjust.Omega} and \code{adjust.n}.
#' @param df [logical] should the degree of freedoms of the Wald statistic be computed using the Satterthwaite correction?
#' Otherwise the degree of freedoms are set to \code{Inf}, i.e. a normal distribution is used instead of a Student's t distribution when computing the p-values.
#' @param adjust.Omega [logical] should the standard errors of the coefficients be corrected for small sample bias?
#' @param adjust.n [logical] should the correction for the degree of freedom be performed?
#' @param tol [numeric >0] the minimum absolute difference between two estimation of the small sample bias.
#' Below this value, the algorithm used to estimate the bias stop.
#' @param n.iter [integer >0] the maximum number of iterations used to estimate the small sample bias of the residual variance-covariance matrix. 
#' @param numeric.derivative [logical] should a numerical derivative be used to compute the first derivative of the information matrix?
#' Otherwise an analytic formula is used.
#' @param trace [logical] should the execution of the function be traced.
#' @param ... [internal] only used by the generic method or by the <- methods.
#'
#' @details The argument \code{value} is equivalent to the argument \code{bias.correct} of the function \code{summary2}.
#' 
#' @concept small sample inference
#' @concept derivative of the score equation
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' dW <- sampleRepeated(10, format = "wide")
#' set.seed(10)
#' dL <- sampleRepeated(10, format = "long")
#' dL$time2 <- paste0("visit",dL$time)
#' 
#' #### linear model ####
#' e.lm <- lm(Y1~X1, data = dW)
#' e2.lm <- sCorrect(e.lm)
#'
#' #### gls model ####
#' e.gls1 <- gls(Y1~X1, data = dW)
#' .getGroups2(e.gls1, data = dW)
#' 
#' e.gls2 <- gls(Y~X1, correlation = corCompSymm(form=~1|id), data = dL)
#' .getGroups2(e.gls2, data = dL)
#'
#' e.gls3 <- gls(Y~X1, correlation = corCompSymm(form=~time|id), data = dL)
#' .getGroups2(e.gls3, data = dL)
#'
#' e.gls4 <- gls(Y~X1, weight = varIdent(form=~1|time2), data = dL)
#' .getGroups2(e.gls4, data = dL)
#'
#' e.gls5 <- gls(Y~X1, weight = varIdent(form=~1|time2),
#'               correlation = corSymm(form=~time|id), data = dL)
#' .getGroups2(e.gls5, data = dL)
#'
#' #### lme model ####
#' e.lme <- lme(Y~X1, random=~1|id, data = dL)
#' .getGroups2(e.lme, data = dL)
#' 
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' sCorrect(e.lvm) <- TRUE ## i.e. bias.correct = TRUE
#' summary2(e.lvm)
#' 
#' @export
`sCorrect` <-
    function(object, param, data,
             first.order, ssc, df,
             derivative, iter.max, tol.max, trace) UseMethod("sCorrect")

## * sCorrect.lm
#' @rdname sCorrect
#' @export
sCorrect.lm <- function(object, param = NULL, data = NULL,
                        first.order = TRUE, ssc = lava.options()$ssc, df = lava.options()$df,
                        derivative = "analytic", iter.max = 100, tol.max = 1e-6,
                        trace = 1){

    ## ** preliminary tests
    if(!inherits(object,"lm") && !inherits(object,"gls") && !inherits(object,"lme") && !inherits(object,"lvmfit")){
        stop("Cannot only convert object of class lm/gls/lme/lvmfit to sCorrect \n")
    }

    ## lvm
    if("multigroupfit" %in% class(object)){
        stop("sCorrect cannot handle multigroup models \n")
    }
    if(inherits(object,"lvmfit") && length(object$model$attributes$ordinal)>0){
        name.t <- names(object$model$attributes$ordinal)
        stop("sCorrect does not handle ordinal variables \n",
             "ordinal variable(s): \"",paste(name.t, collapse = "\" \""),"\"\n")
    }
    
    if(inherits(object,"lvmfit") && length(object$model$attributes$transform)>0){
        name.t <- names(object$model$attributes$transform)
        stop("sCorrect does not handle transformed variables \n",
             "transformed variable(s): \"",paste(name.t, collapse = "\" \""),"\"\n")
    }

    ## arguments
    if(!is.null(ssc) && ssc %in% c("residuals","Cox") == FALSE){
        stop("Argument \'ssc\' should either be NULL, \"residual\", or \"Cox\" \n")
    }
    if(!is.null(df) && df %in% c("Satterthwaite") == FALSE){
        stop("Argument \'df\' should either be NULL or \"Satterthwaite\" \n")
    }
    
    ## ** initialize object
    if(trace>0){cat("Initialize moments \n")}
    object <- .init_sCorrect(object, param = param, data = data,
                             initialize = TRUE, 
                             first.order = first.order, dVcov = !is.null(df), dVcov.robust = !is.null(df),
                             derivative = derivative)

    ## ** bias correction
    if(!is.null(ssc)){

        ## *** initialize bias correction
        if(ssc=="Cox"){
            iParam <- .init_sscCoxSnell(object, param = iParam)
        }else if(ssc=="residuals"){
            iParam <- .init_sscResiduals(object, param = iParam)
        }
        
        ## *** perform bias correction
        if(trace>0){cat("Perform bias correction \n")}
        iIter <- 0
        iTol <- Inf
        iParam <- object$sCorrect$coef
        iiParam <- iParam
        
        while(iIter < iter.max & iTol > tol.max){
            if(trace>0){cat("*")}

            ## bias correction
            if(ssc=="Cox"){
                iParam <- .sscCoxSnell(object, param = iParam)
            }else if(ssc=="residuals"){
                iParam <- .sscResiduals(object, param = iParam)
            }
            
            ## update moments
            object <- .init_sCorrect(object, param = iParam,
                                     initialize = FALSE,
                                     first.order = first.order, dVcov = (ssc=="Cox"), dVcov.robust = (ssc=="Cox"),
                                     derivative = derivative)

            ## cv criteria
            iIter <- iIter + 1
            iTol <- max(abs(iParam-iiParam))
            iiParam <- iParam
        }

        ## *** assess convergence
        if(iTol > tol.max){
            cv <- FALSE
            if(trace > 0){
                warning("small sample correction did not reach convergence after ",iIter," iterations \n")
            }
        }else{
            cv <- TRUE
            if(trace > 0){
                cat("\n")
            }
        }
        
        ## *** update score, information matrix, leverage, ..., with the bias corrected parameters
        object <- .init_sCorrect(object, param = iParam,
                                 initialize = FALSE,
                                 first.order = first.order, dVcov = !is.null(df), dVcov.robust = !is.null(df),
                                 derivative = derivative)
    }
    
    ## ** export    
    object$sCorrect$ssc <- ssc
    if(!is.null(ssc)){
        object$sCorrect$opt.scc <- list(cv = cv, iter = iIter, tol = iTol, iter.max = iter.max, tol.max = tol.max)
    }
    object$sCorrect$df <- df
    class(object) <- append("sCorrect",class(object))
    return(object)    
}

## * sCorrect.gls
#' @rdname sCorrect
#' @export
sCorrect.gls <- sCorrect.lm

## * sCorrect.lme
#' @rdname sCorrect
#' @export
sCorrect.lme <- sCorrect.lm

## * sCorrect.lvmfit
#' @rdname sCorrect
#' @export
sCorrect.lvmfit <- sCorrect.lm

## * sCorrect.sCorrect
#' @rdname sCorrect
#' @export
sCorrect.sCorrect <- function(object, param = NULL, data = NULL,
                              first.order = TRUE, ssc = lava.options()$ssc, df = lava.options()$df,
                              derivative = "analytic", iter.max = 100, tol.max = 1e-6,
                              trace = 1){
    class(object) <- setdiff(class(object),"sCorrect")
    return(sCorrect(object, param = param, data = data,
                    first.order = first.order, ssc = ssc, df = df,
                    derivative = derivative, iter.max = iter.max, tol.max = tol.max,
                    trace = trace))    
}

## * .init_sCorrect
.init_sCorrect <- function(object, param, data,
                           initialize,
                           first.order = first.order, dVcov, dVcov.robust,
                           derivative){
    if(lava.options()$debug){cat(".init_sCorrect \n")}

    derivative <- match.arg(derivative, choices = c("analytic","numeric"))

    ## ** initalize sCorrect object with model parameters, data, ...
    object$sCorrect <- conditionalMoment(object,
                                         initialize = initialize, first.order = first.order,
                                         second.order = dVcov || dVcov.robust,
                                         usefit = TRUE,
                                         param = param,
                                         data = data)

    if(initialize){ ## leverage
        object$sCorrect$leverage <- matrix(0, nrow = object$sCorrect$cluster$n.cluster, ncol = length(object$sCorrect$endogenous),
                                           dimnames = list(NULL, object$sCorrect$endogenous))
        object$sCorrect$leverage[is.na(object$sCorrect$moment$residuals)] <- NA
    }
    
    ## ** score
    if(first.order){
        object$sCorrect$score <- .score2(dmu = object$sCorrect$dmoment$dmu,
                                         dOmega = object$sCorrect$dmoment$dOmega,                    
                                         epsilon = object$sCorrect$moment$residuals,
                                         OmegaM1 = object$sCorrect$moment$OmegaM1.missing.pattern,
                                         missing.pattern = object$sCorrect$missing$pattern,
                                         unique.pattern = object$sCorrect$missing$unique.pattern,
                                         name.pattern = object$sCorrect$missing$name.pattern,
                                         name.param = object$sCorrect$skeleton$Uparam,
                                         name.meanparam = object$sCorrect$skeleton$Uparam.mean,
                                         name.varparam = object$sCorrect$skeleton$Uparam.var,
                                         n.cluster = object$sCorrect$cluster$n.cluster)
    }
    
    ## ** information matrix
    if(first.order){
        object$sCorrect$information <- .information2(dmu = object$sCorrect$dmoment$dmu,
                                                     dOmega = object$sCorrect$dmoment$dOmega,
                                                     OmegaM1 = object$sCorrect$moment$OmegaM1.missing.pattern,
                                                     missing.pattern = object$sCorrect$missing$pattern,
                                                     unique.pattern = object$sCorrect$missing$unique.pattern,
                                                     name.pattern = object$sCorrect$missing$name.pattern,
                                                     grid.mean = object$sCorrect$skeleton$grid.dmoment$mean, 
                                                     grid.var = object$sCorrect$skeleton$grid.dmoment$var, 
                                                     name.param = object$sCorrect$skeleton$Uparam,
                                                     leverage = object$sCorrect$leverage,
                                                     n.cluster = object$sCorrect$cluster$n.cluster)
        object$sCorrect$vcov.param  <- .info2vcov(object$sCorrect$information, attr.info = FALSE)
    }
    
    ## ** hessian (analytic)
    if((derivative == "analytic") && dVcov.robust){
        object$sCorrect$hessian <- .hessian2(dmu = object$sCorrect$dmoment$dmu,
                                             dOmega = object$sCorrect$dmoment$dOmega,
                                             d2mu = object$sCorrect$d2moment$d2mu,
                                             d2Omega = object$sCorrect$d2moment$d2Omega,
                                             epsilon = object$sCorrect$moment$residuals,                                     
                                             OmegaM1 = object$sCorrect$moment$OmegaM1.missing.pattern,
                                             missing.pattern = object$sCorrect$missing$pattern,
                                             unique.pattern = object$sCorrect$missing$unique.pattern,
                                             name.pattern = object$sCorrect$missing$name.pattern,
                                             grid.mean = object$sCorrect$skeleton$grid.dmoment$mean, 
                                             grid.var = object$sCorrect$skeleton$grid.dmoment$var, 
                                             grid.hybrid = object$sCorrect$skeleton$grid.dmoment$hybrid, 
                                             name.param = object$sCorrect$skeleton$Uparam,
                                             leverage = object$sCorrect$leverage,
                                             n.cluster = object$sCorrect$cluster$n.cluster)
    }

    ## ** dVcov.param (model based variance, analytic)
    if((derivative == "analytic") && (dVcov || dVcov.robust)){
        object$sCorrect$dInformation <- .dInformation2(dmu = object$sCorrect$dmoment$dmu,
                                                       dOmega = object$sCorrect$dmoment$dOmega,
                                                       d2mu = object$sCorrect$d2moment$d2mu,
                                                       d2Omega = object$sCorrect$d2moment$d2Omega,
                                                       OmegaM1 = object$sCorrect$moment$OmegaM1.missing.pattern,
                                                       missing.pattern = object$sCorrect$missing$pattern,
                                                       unique.pattern = object$sCorrect$missing$unique.pattern,
                                                       name.pattern = object$sCorrect$missing$name.pattern,
                                                       grid.dInformation = object$sCorrect$skeleton$grid.dInformation,
                                                       name.param = object$sCorrect$skeleton$Uparam,
                                                       name.param.dInformation = object$sCorrect$skeleton$name.param.dInformation,
                                                       leverage = object$sCorrect$leverage,
                                                       n.cluster = object$sCorrect$cluster$n.cluster)

        ## delta method 
        object$sCorrect$dVcov.param <- array(NA, dim = dim(object$sCorrect$dInformation), dimnames = dimnames(object$sCorrect$dInformation))
        for(iP in object$sCorrect$skeleton$name.param.dInformation){
            object$sCorrect$dVcov.param[,,iP] <- - object$sCorrect$vcov.param %*% object$sCorrect$dInformation[,,iP] %*% object$sCorrect$vcov.param
        }
    }

    ## ** dRvcov.param  (robust variance, analytic)
    if((derivative == "analytic") && dVcov.robust){
        name.param <- object$sCorrect$skeleton$Uparam
        name.param.dInformation <- object$sCorrect$skeleton$name.param.dInformation
        n.param <- length(name.param)

        object$sCorrect$dRvcov.param <- array(NA, dim = c(n.param,n.param,n.param), dimnames = list(name.param,name.param,name.param))
        score <- object$sCorrect$score
        score2 <- crossprod(score)
        hessian <- object$sCorrect$hessian
        vcov.param <- object$sCorrect$vcov.param
        dVcov.param <- object$sCorrect$dVcov.param
        
        for(iP in 1:n.param){ ## iP <- 1
            if(name.param[iP] %in% name.param.dInformation){
                term1 <- dVcov.param[,,name.param[iP]] %*% score2 %*% vcov.param
            }else{
                term1 <- matrix(0, nrow = n.param, ncol = n.param)
            }
            term2 <- vcov.param %*% hessian[iP,,] %*% score %*% vcov.param
            object$sCorrect$dRvcov.param[,,iP] <- term1 + t(term1) + term2 + t(term2)
        }
    }

    ## ** dVcov.param and dRvcov.param (numeric derivatives)
    if(derivative == "numeric" && (dVcov || dVcov.robust)){
        test.package <- try(requireNamespace("numDeriv"), silent = TRUE)
        if(inherits(test.package,"try-error")){
            stop("There is no package \'numDeriv\' \n",
                 "This package is necessary when argument \'numeric.derivative\' is TRUE \n")
        }

        ## range(object$sCorrect$score - .warper.numDev(value = object$sCorrect$param, object = object, type = "score"))
        ## range(object$sCorrect$hessian - .warper.numDev(value = object$sCorrect$param, object = object, type = "hessian"))
        ## range(object$sCorrect$vcov.param - .warper.numDev(value = object$sCorrect$param, object = object, type = "vcov.model"))

        param <- object$sCorrect$param
        name.param <- names(object$sCorrect$param)
        n.param <- length(param)
        n.cluster <- object$sCorrect$cluster$n.cluster

        ## *** hessian
        num.hessian <- numDeriv::jacobian(.warper.numDev, x = param, object = object, type = "score", method = "Richardson")

        object$sCorrect$hessian <- aperm(array(num.hessian, dim = c(n.cluster,n.param,n.param),
                                               dimnames = list(NULL, name.param, name.param)), perm = 3:1)

        ## ## *** dInformation
        ## num.information <- numDeriv::jacobian(.warper.numDev, x = param, object = object, type = "information", method = "Richardson")

        ## object$sCorrect$dInformation <- array(num.information, dim = c(n.param,n.param,n.param),
        ##                                       dimnames = list(name.param, name.param, name.param))

        ## *** dVcov.param
        num.dVcov.param <- numDeriv::jacobian(.warper.numDev, x = param, object = object, type = "vcov.model", method = "Richardson")

        object$sCorrect$dVcov.param <- array(num.dVcov.param, dim = c(n.param,n.param,n.param),
                                             dimnames = list(name.param, name.param, name.param))


        ## *** dRvcov.param
        num.dRvcov.param <- numDeriv::jacobian(.warper.numDev, x = param, object = object, type = "vcov.robust", method = "Richardson")

        object$sCorrect$dRvcov.param <- array(num.dRvcov.param, dim = c(n.param,n.param,n.param),
                                              dimnames = list(name.param, name.param, name.param))
    }        
        
    ## ** output
    return(object)
}

## * .wraper.numDev
.warper.numDev <- function(value, object, type){ # x <- p.obj

    type <- match.arg(type, c("score","hessian","information","vcov.model","vcov.robust"))
            
    ## update param
    pp <- object$sCorrect$coef
    pp[names(value)] <- value

    ## update moments and their derivatives
    cM <- conditionalMoment(object, param = pp,
                            initialize = FALSE, first.order = TRUE, second.order = (type == "hessian"), usefit = TRUE)

    ## compute information matrix
    if(type %in% c("score","vcov.robust")){
        ss <- .score2(dmu = cM$dmoment$dmu,
                      dOmega = cM$dmoment$dOmega,
                      epsilon = cM$moment$residuals,
                      OmegaM1 = cM$moment$OmegaM1.missing.pattern,
                      missing.pattern = object$sCorrect$missing$pattern,
                      unique.pattern = object$sCorrect$missing$unique.pattern,
                      name.pattern = object$sCorrect$missing$name.pattern,
                      name.param = object$sCorrect$skeleton$Uparam,
                      name.meanparam = object$sCorrect$skeleton$Uparam.mean,
                      name.varparam = object$sCorrect$skeleton$Uparam.var,
                      n.cluster = object$sCorrect$cluster$n.cluster)
    }
    if(type=="hessian"){
        hh <- .hessian2(dmu = cM$dmoment$dmu,
                        dOmega = cM$dmoment$dOmega,
                        d2mu = cM$d2moment$d2mu,
                        d2Omega = cM$d2moment$d2Omega,
                        epsilon = cM$moment$residuals,
                        OmegaM1 = cM$moment$OmegaM1.missing.pattern,
                        missing.pattern = object$sCorrect$missing$pattern,
                        unique.pattern = object$sCorrect$missing$unique.pattern,
                        name.pattern = object$sCorrect$missing$name.pattern,
                        grid.mean = object$sCorrect$skeleton$grid.dmoment$mean,
                        grid.var = object$sCorrect$skeleton$grid.dmoment$var,
                        grid.hybrid = object$sCorrect$skeleton$grid.dmoment$hybrid,
                        name.param = object$sCorrect$skeleton$Uparam,
                        leverage = object$sCorrect$leverage,
                        n.cluster = object$sCorrect$cluster$n.cluster)
    }
    if(type %in% c("information","vcov.model","vcov.robust")){
        ii <- .information2(dmu = cM$dmoment$dmu,
                            dOmega = cM$dmoment$dOmega,
                            OmegaM1 = cM$moment$OmegaM1.missing.pattern,
                            missing.pattern = object$sCorrect$missing$pattern,
                            unique.pattern = object$sCorrect$missing$unique.pattern,
                            name.pattern = object$sCorrect$missing$name.pattern,
                            grid.mean = object$sCorrect$skeleton$grid.dmoment$mean,
                            grid.var = object$sCorrect$skeleton$grid.dmoment$var,
                            name.param = object$sCorrect$skeleton$Uparam,
                            leverage = object$sCorrect$leverage,
                            n.cluster = object$sCorrect$cluster$n.cluster)
    }

    return(switch(type,
                  "score" = ss,
                  "hessian" = hh,
                  "information" = ii,
                  "vcov.model" = solve(ii),
                  "vcov.robust" = crossprod(ss %*% solve(ii)))
           )
}
     
##----------------------------------------------------------------------
### sCorrect.R ends here










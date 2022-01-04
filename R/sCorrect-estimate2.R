### sCorrect-estimate2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2018 (14:29) 
## Version: 
## Last-Updated: Jan  4 2022 (16:50) 
##           By: Brice Ozenne
##     Update #: 2029
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - estimate2
#' @title  Satterthwaite Correction and Small Sample Correction
#' @description Correct the bias of the ML estimate of the variance and compute the first derivative of the information matrix.
#' @name estimate2
#'
#' @param object a \code{lvm} object.
#' @param param [numeric vector, optional] the values of the parameters at which to perform the correction.
#' @param data [data.frame, optional] the dataset relative to which the correction should be performed.
#' @param ssc [character] type of small sample correction: \code{"cox"} or \code{"residuals"}, or \code{"none"}.
#' @param df [character] method used to estimate the degree of freedoms of the Wald statistic: \code{"satterthwaite"}, or \code{"none"}
#' Otherwise the degree of freedoms are set to \code{Inf}, i.e. a normal distribution is used instead of a Student's t distribution when computing the p-values.
#' @param tol.max [numeric >0] the largest acceptable absolute difference between two succesive estimates of the bias correction.
#' @param iter.max [integer >0] the maximum number of iterations used to estimate the bias correction.
#' @param derivative [character] should the first derivative of the information matrix be computed using a formula (\code{"analytic"}) or numerical derivative (\code{"numeric"})?
#' @param hessian [logical] should the hessian be stored? Can be \code{NULL} to indicate only if computed during the small sample correction.
#' @param trace [logical] should the execution of the function be traced.
#' @param ...  arguments passed to \code{lava::estimate} when using a \code{lvm} object.
#'
#' @details The argument \code{value} is equivalent to the argument \code{bias.correct} of the function \code{summary2}.
#' 
#' @concept estimator
#' @keywords smallSampleCorrection
#' 
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' dW <- sampleRepeated(10, format = "wide")
#' 
#' #### latent variable model ####
#' m.lvm <- lvm(Y1~X1+X2+Z1)
#'
#' e2.lvm <- estimate2(m.lvm, data = dW)
#' summary2(e2.lvm)
#' 
#' @export
`estimate2` <-
    function(object, param, data,
             ssc, df,
             derivative, hessian, iter.max, tol.max, trace , ...) UseMethod("estimate2")


## * estimate2.lvm
#' @rdname estimate2
estimate2.lvm <- function(object, data = NULL,
                          ssc = lava.options()$ssc, df = lava.options()$df,
                          derivative = "analytic", hessian = NULL, iter.max = 100, tol.max = 1e-6,
                          trace = 0, ...){

    out <- lava::estimate(x = object, data = data, ...)
    return(estimate2(out, param = NULL, data = NULL,
                     ssc = ssc, df = df,
                     derivative = derivative, hessian = hessian, iter.max = iter.max, tol.max = tol.max, trace = trace))

}

## * estimate2.lvmfit
#' @rdname estimate2
estimate2.lvmfit <- function(object, param = NULL, data = NULL,
                             ssc = lava.options()$ssc, df = lava.options()$df,
                             derivative = "analytic", hessian = NULL, iter.max = 100, tol.max = 1e-6,
                             trace = 0, ...){

    ## ** preliminary tests
    if(length(list(...))>0){
        stop("Argument \"",paste(names(list(...)), collapse = "\" \""),"\" not used. \n")
    }

    if("multigroupfit" %in% class(object)){
        stop("estimate2 cannot handle multigroup models \n")
    }

    if(inherits(object,"lvmfit") && length(object$model$attributes$ordinal)>0){
        name.t <- names(object$model$attributes$ordinal)
        stop("estimate2 does not handle ordinal variables \n",
             "ordinal variable(s): \"",paste(name.t, collapse = "\" \""),"\"\n")
    }
    
    if(inherits(object,"lvmfit") && length(object$model$attributes$transform)>0){
        name.t <- names(object$model$attributes$transform)
        stop("estimate2 does not handle transformed variables \n",
             "transformed variable(s): \"",paste(name.t, collapse = "\" \""),"\"\n")
    }

    ## arguments
    ssc <- match.arg(tolower(ssc), c("none","residuals","cox"))
    df <- match.arg(tolower(df), c("none","satterthwaite"))

    if(df %in% "satterthwaite" || ssc %in% "cox"){
        second.order <- TRUE
    }else if(ssc %in% "residuals"){
        second.order <- FALSE
    }else{
        second.order <- FALSE
    }
    
    ## ** initialize object
    if(trace>0){cat("Initialization:")}
    object$sCorrect <- moments2(object, data = data, param = param, Psi = NULL,
                                initialize = TRUE, usefit = TRUE,
                                score = TRUE, information = TRUE, hessian = hessian, vcov = TRUE,
                                dVcov = (ssc == "cox")  || (ssc == "none" && df == "satterthwaite"), dVcov.robust = FALSE,
                                residuals = TRUE, leverage = FALSE, derivative = derivative)  ## setting leverage to FALSE is like initialization to 0
    if(trace>0){cat(" done \n")}

    ## ** bias correction    
    if(ssc != "none"){

        ## *** initialize bias correction
        if(trace>0){cat("Initialize bias correction \n")}
        if(ssc=="Cox"){
            object.ssc <- list(type = "Cox",
                               param0 = object$sCorrect$param,
                               Omega0 = object$sCorrect$moment$Omega)
        }else if(ssc=="residuals"){
            object.ssc <- .init_sscResiduals(object)
        }
        
        ## *** perform bias correction
        if(trace>0){cat("Perform bias correction \n")}
        iCV <- FALSE
        iIter <- 0
        iTol <- Inf
        iiParam <- object$sCorrect$param
        ## cat(iTol," (",iiParam,") \n")
            
        while(iCV == FALSE && iIter < iter.max){
            if(trace>0){cat("*")}

            ## bias correction
            if(ssc=="Cox"){
                iParam <- .sscCoxSnell(object, ssc = object.ssc)
                object.ssc$JJK <- attr(iParam,"JJK")
                object.ssc$lm <- attr(iParam,"lm")
                object.ssc$Psi <- object$sCorrect$moment$Omega - object.ssc$Omega0
            }else if(ssc=="residuals"){
                iParam <- .sscResiduals(object, ssc = object.ssc)
                object.ssc$Psi <- attr(iParam,"Psi")
                object$sCorrect$leverage <- attr(iParam,"leverage")
                object$sCorrect$residuals <- attr(iParam,"residuals")
            }
            ## update moments
            ## object.ssc$Omega0 + object.ssc$Psi - object$sCorrect$moment$Omega
            object$sCorrect <- moments2(object, param = iParam, Psi = object.ssc$Psi, 
                                        initialize = FALSE, usefit = TRUE,
                                        score = TRUE, information = TRUE, hessian = NULL, vcov = TRUE,
                                        dVcov = (ssc == "cox"), dVcov.robust = FALSE,
                                        residuals = TRUE, leverage = TRUE, derivative = derivative)

            ## cv criteria
            iIter <- iIter + 1
            iTol <- max(abs(iParam-iiParam))
            ## cat(iTol," (",iParam,") \n")
            iiParam <- iParam
            iCV <- iTol <= tol.max
        }

        ## *** assess convergence
        object.ssc$cv <- iCV
        object.ssc$iter <- iIter
        object.ssc$tol <- iTol
        object.ssc$iter.max <- iter.max
        object.ssc$tol.max <- tol.max
        
        if(iCV == FALSE  && trace > 0){
            warning("small sample correction did not reach convergence after ",iIter," iterations \n")
        }
        if(trace > 0){
            cat("\n")
        }

        ## *** update score, information matrix, leverage, ..., with the bias corrected parameters
        if(trace>0){cat("Update moments/df \n")}
        object$sCorrect <- moments2(object, data = data, param = iParam, Psi = object.ssc$Psi,
                                    initialize = FALSE, usefit = TRUE,
                                    score = TRUE, information = TRUE, hessian = hessian, vcov = TRUE,
                                    dVcov = (df == "satterthwaite"), dVcov.robust = FALSE,
                                    residuals = TRUE, leverage = TRUE, derivative = derivative)
        object$sCorrect$ssc <- object.ssc
    }else{
        object$sCorrect$ssc$type <- NA
    }

    ## ** degrees of freedom    
    object$sCorrect$df <- df ## degrees of freedom are computed later (by compare2)
    
    ## ** restaure original param order
    name.param <- object$sCorrect$name.param
    if(!is.null(name.param)){
        object$sCorrect$param <- object$sCorrect$param[name.param]
        names(object$sCorrect$param) <- names(name.param)

        if(!is.null(object$sCorrect$score)){
            object$sCorrect$score <- object$sCorrect$score[,name.param,drop=FALSE]
            colnames(object$sCorrect$score) <- names(name.param)
        }
        if(!is.null(object$sCorrect$information)){
            object$sCorrect$information <- object$sCorrect$information[name.param,name.param,drop=FALSE]
            dimnames(object$sCorrect$information) <- list(names(name.param),names(name.param))
        }
        if(!is.null(object$sCorrect$vcov.param)){
            object$sCorrect$vcov.param <- object$sCorrect$vcov.param[name.param,name.param,drop=FALSE]
            dimnames(object$sCorrect$vcov.param) <- list(names(name.param),names(name.param))
        }
        if(!is.null(object$sCorrect$hessian)){
            object$sCorrect$hessian <- object$sCorrect$hessian[name.param,name.param,,drop=FALSE]
            dimnames(object$sCorrect$hessian) <- list(names(name.param),names(name.param),NULL)
        }
        if(!is.null(object$sCorrect$dInformation)){
            object$sCorrect$dInformation <- object$sCorrect$dInformation[name.param,name.param,name.param,drop=FALSE]
            dimnames(object$sCorrect$dInformation) <- list(names(name.param),names(name.param),names(name.param))
        }
        if(!is.null(object$sCorrect$dVcov.param)){
            object$sCorrect$dVcov.param <- object$sCorrect$dVcov.param[name.param,name.param,name.param,drop=FALSE]
            dimnames(object$sCorrect$dVcov.param) <- list(names(name.param),names(name.param),names(name.param))
        }
        if(!is.null(object$sCorrect$dRvcov.param)){
            object$sCorrect$dRvcov.param <- object$sCorrect$dRvcov.param[name.param,name.param,name.param,drop=FALSE]
            dimnames(object$sCorrect$dRvcov.param) <- list(names(name.param),names(name.param),names(name.param))
        }
    }
    
    ## ** export
    class(object) <- append("lvmfit2",class(object))
    return(object)    
}

## * estimate2.list
#' @rdname estimate2
estimate2.list <- function(object, ...){
    object.class <- class(object)
    object <- lapply(object, sCorrect, ...)
    class(object) <- object.class
    return(object)
}

## * estimate2.mmm
#' @rdname estimate2
estimate2.mmm <- estimate2.list

##----------------------------------------------------------------------
### sCorrect.R ends here










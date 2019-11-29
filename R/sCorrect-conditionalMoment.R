### conditionalMoment.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: nov 29 2019 (18:17) 
##           By: Brice Ozenne
##     Update #: 1224
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * conditionalMoment - documentation
#' @title Prepare the Computation of score2
#' @description Compute the conditional mean and variance,
#' and their first and second derivative regarding the model parameters.
#' @name conditionalMoment
#' 
#' @param object,x a latent variable model.
#' @param data [data.frame] data set.
#' @param param [numeric vector] the fitted coefficients.
#' @param cluster [list] information about the cluster and endogenous variables (output of getGroups2).
#' @param first.order [logical] should the terms relative to the first derivative of the conditional moments be considered?
#' @param second.order [logical] should the terms relative to the second derivative of the conditional moments be considered?
#' @param usefit [logical] If TRUE the coefficients estimated by the model are used to pre-compute quantities. 
#' @param ... [internal] only used by the generic method or by the <- methods.
#' 
#' @details For lvmfit objects, there are two levels of pre-computation:
#' \itemize{
#' \item a basic one that do no involve the model coefficient (\code{conditionalMoment.lvm}).
#' \item an advanced one that require the model coefficients (\code{conditionalMoment.lvmfit}). 
#' }
#' 
#' @examples
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' d <- lava::sim(m,1e2)
#' e <- estimate(m, d)
#'
#' ## basic pre-computation
#' res1 <- conditionalMoment(e, data = d,
#'                          first.order = FALSE, second.order = FALSE,
#'                          name.endogenous = endogenous(e),
#'                          name.latent = latent(e), usefit = FALSE)
#' res1$skeleton$Sigma
#' 
#' ## full pre-computation
#' res2 <- conditionalMoment(e, param = coef(e), data = d,
#'                          first.order = FALSE, second.order = FALSE,
#'                          name.endogenous = endogenous(e),
#'                          name.latent = latent(e), usefit = TRUE
#' )
#' res2$value$Sigma
#'
#' @concept small sample inference
#' @concept derivative of the score equation
#' 
#' @keywords internal
#' @export
`conditionalMoment` <-
  function(object, ...) UseMethod("conditionalMoment")

## * conditionalMoment.lm
#' @rdname conditionalMoment
#' @export
conditionalMoment.lm <- function(object,
                                 initialize, first.order, second.order, usefit,
                                 param = NULL,
                                 X = NULL,
                                 cluster = NULL,
                                 endogenous = NULL,
                                 latent = NULL,
                                 ...){
    if(TRUE){cat("conditionalMoment \n")}
    ## ** get information from object
    if(is.null(endogenous)){
        endogenous <- lava::endogenous(object, format = "wide")
    }
    if(is.null(latent)){
        latent <- lava::latent(object)
    }
    if(is.null(X)){
        data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                            envir = parent.env(environment()), rm.na = TRUE)
        X <- .getDesign(object, data = data, add.Y = TRUE)
    }
    if(is.null(cluster)){
        data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                            envir = parent.env(environment()), rm.na = TRUE)
        cluster <- .getGroups2(object, data = data, endogenous = endogenous)
    }

    if(any(colnames(X) %in% c("XXvalueXX","XXendogenousXX","XXclusterXX"))){
        stop("\"XXvalueXX\", \"XXendogenousXX\", and \"XXclusterXX\" should not correspond to variable names \n",
             "It is used internally for data manipulation \n")
    }

    if(NROW(X)==length(cluster$name.cluster)){## convert to long format
        X.long <- melt(as.data.frame(X), id.vars = setdiff(colnames(X),endogenous),
                       measure.vars = endogenous,
                       variable.name = "XXendogenousXX",
                       value.name = "XXvalueXX")
    }else{
        X.long <- X
        names(X.long)[names(X.long) == lava::endogenous(object, format = "long")] <- "XXvalueXX"
        X.long$XXendogenousXX <- NA
        for(iC in 1:cluster$n.cluster){ ## iC <- 10
            X.long$XXendogenousXX[cluster$index.cluster == iC] <- na.omit(cluster$index.Omega[[iC]])
        }
    }
    X.long <- cbind(X.long, XXclusterXX = cluster$index.cluster)
    X.long <- X.long[order(X.long$XXclusterXX,X.long$XXendogenousXX),]

    ## ** sanity checks
    if(first.order == FALSE && second.order == TRUE){
        stop("Cannot pre-compute second order derivatives without the first order derivatives. \n")
    }
    if(initialize == FALSE && is.null(object$sCorrect$conditionalMoment)){
        stop("Initialization of the conditional moments missing. \n",
             "Consider setting the argument \'initialize\' to TRUE \n")
    }

    ## ** initialize
    if(initialize){
        object$sCorrect$conditionalMoment <- list(moment = NULL,
                                                  dMoment = NULL,
                                                  d2Moment = NULL)

        ## *** initialize conditional moments
        object$sCorrect$conditionalMoment$moment <- initSkeleton(object, X = X.long,
                                                                 endogenous = endogenous, latent = latent,
                                                                 as.lava = TRUE)

        ## *** initialize partial derivatives of the conditional moments
        if(first.order){
            object$sCorrect$conditionalMoment$dmoment <- initSkeletonDtheta(object, data = X,
                                                                            df.param.all = object$sCorrect$conditionalMoment$moment$df.param,
                                                                            param2originalLink = object$sCorrect$conditionalMoment$moment$param2originalLink,
                                                                            endogenous = endogenous, 
                                                                            latent = latent)
        }
    
        ## *** initialize second order partial derivatives of the conditional moments
        if(second.order){
            object$sCorrect$conditionalMoment$d2moment <- initSkeletonDtheta2(object, data = X,
                                                                              df.param.all = object$sCorrect$conditionalMoment$moment$df.param,
                                                                              param2originalLink = object$sCorrect$conditionalMoment$moment$param2originalLink,
                                                                              latent = latent)
        }

        ## *** param with non-zero third derivative
        type.3deriv <- c("alpha","Gamma","Lambda","B","Psi_var","Sigma_var","Psi_cov","Sigma_cov")
        index.keep <- intersect(which(!is.na(object$sCorrect$conditionalMoment$df.param$lava)),
                                which(object$sCorrect$conditionalMoment$df.param$detail %in% type.3deriv)
                                )    
        object$sCorrect$conditionalMoment$name.3deriv <- object$sCorrect$conditionalMoment$df.param[index.keep, "originalLink"]
    }

    ## ** update according to the value of the model coefficients
    if(usefit){

        if(is.null(param)){
            param <- coef2(object)
        }

        ## *** conditional moments
        object$sCorrect$conditionalMoment$value <- skeleton(object, data = X, param = param,
                                                            endogenous = endogenous, 
                                                            latent = latent)

        if(object$sCorrect$conditionalMoment$skeleton$toUpdate["param"]){
            object$sCorrect$conditionalMoment$param <- coef(object)
        }
        if(object$sCorrect$conditionalMoment$skeleton$toUpdate["mu"]){            
            if(n.latent==0){
                object$sCorrect$conditionalMoment$mu <- object$sCorrect$conditionalMoment$value$nu.XK
            }else{
                object$sCorrect$conditionalMoment$mu <- object$sCorrect$conditionalMoment$value$nu.XK + object$sCorrect$conditionalMoment$value$alpha.XGamma.iIB %*% object$sCorrect$conditionalMoment$value$Lambda
            }            
        }
        if(object$conditionalMoment$skeleton$toUpdate["Omega"]){
            object$conditionalMoment$Omega <- getVarCov2(object)
        }
        
        ## *** first order derivatives
        if(first.order){            
            out <- skeletonDtheta(object,
                                  endogenous = endogenous, 
                                  latent = latent)
            object$conditionalMoment$dmu <- out$dmu
            object$conditionalMoment$dOmega <- out$dOmega            
        }

        ## *** second order derivatives
        if(second.order){
            out2 <- skeletonDtheta2(object)
            object$conditionalMoment$d2mu <- out2$d2mu
            object$conditionalMoment$d2Omega <- out2$d2Omega
        }
       
    }
     
    ## ** export
    return(object$conditionalMoment)
}

## * conditionalMoment.gls
#' @rdname conditionalMoment
#' @export
conditionalMoment.gls <- conditionalMoment.lm

## * conditionalMoment.lme
#' @rdname conditionalMoment
#' @export
conditionalMoment.lme <- conditionalMoment.lm

## * conditionalMoment.lvmfit
#' @rdname conditionalMoment
#' @export
conditionalMoment.lvmfit <- conditionalMoment.lm

## * conditionalMoment.sCorrect
#' @rdname conditionalMoment
#' @export
conditionalMoment.sCorrect <- conditionalMoment.lm

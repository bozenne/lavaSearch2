### sCorrect_init.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  2 2019 (11:52) 
## Version: 
## Last-Updated: nov 28 2019 (11:59) 
##           By: Brice Ozenne
##     Update #: 169
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .convert_sCorrect (generate object)
## ** .convert_sCorrect<-
`.convert_sCorrect<-` <-
    function(object, value) UseMethod("convert_sCorrect<-")

## ** .convert_sCorrect<-.default
`convert_sCorrect<-.default` <- function(object, ..., value){

    ## *** preliminary tests
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

    ## *** initialize object
    if(is.null(object$sCorrect) || (value == TRUE)){
        object$sCorrect <- .init_sCorrect(object, ...)
    }else{
        stop("Existing small sample correction \n",
             "Set argument \'value\' to TRUE to overwrite it \n")
    }

    ## *** output
    class(object) <- append("sCorrect",class(object))
    return(object)
}

## ** .init_sCorrect
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
#' .init_sCorrect(e.lm)
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

.init_sCorrect <- function(object, ...){
    if(TRUE){cat(".init_sCorrect \n")}

    out <- list()

    ## ** type of variables
    out$endogenous <- endogenous(object)
    out$exogenous <- exogenous(object)
    out$latent <- latent(object)

    ## ** model coefficients
    out$coef <- coef2(object)

    ## ** data
    out$data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                            envir = parent.env(environment()), rm.na = TRUE)

    ## ** design matrix
    out$X <- .getDesign(object, data = out$data, add.Y = TRUE)

    ## ** clusters
    out$cluster <- .getGroups2(object, data = out$data, endogenous = out$endogenous)

    ## ** sample size
    out$n.corrected <- out$cluster$n.cluster

    ## ** leverage
    out$leverage <- matrix(NA, nrow = out$cluster$n.cluster, ncol = length(out$endogenous),
                           dimnames = list(NULL, out$endogenous))
    if(out$cluster$missing==FALSE){
        out$leverage[] <- 0
    }else{
        for(iC in 1:out$cluster$n.cluster){ # iC <- 1
            iIndex <- which(out$cluster$index.cluster==iC)
            out$leverage[iC,iIndex[out$cluster$index.Omega[[iC]]]] <- 0
        }
    }

    ## moments
    out$moment <- conditionalMoment(object,
                                    param = out$coef,
                                    X = out$X,
                                    cluster = out$cluster,
                                    endogenous = out$endogenous,
                                    latent = out$latent,
                                    first.order = TRUE, second.order = TRUE, usefit = FALSE)
    browser()
    
    ## residuals
    out$residuals <- residuals2(object)

    ## out$score <- score2(object)
    out$vcov.param <- vcov2(object)
    ## out$vcov.residual <- getVarCov2(object)
    out <- .information2(dmu,
                         dOmega,
                         Omega,
                         grid.meanparam,
                         n.grid.meanparam,
                         grid.varparam,
                         n.grid.varparam,
                         name.param,
                         n.param)

    
    
    ## *** exclude clusters containing only missing values

    ## *** output
    return(out)
}

## * information2
#' @title  Extract The Full Information Matrix After Small Sample Correction. 
#' @description  Extract the full information matrix from  from \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} objects.
#' Similar to \code{lava::information} but with small sample correction.
#' @name information2
#'
#' @param object a linear model or a latent variable model
#' @param ... arguments to be passed to \code{vcov2}.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @return A matrix.
#' 
#' @examples
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#'
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- sequence.lvm(0)
#' set.seed(10)
#' d <- lava::sim(m,n)
#'
#' ## linear model
#' e.lm <- lm(formula.lvm,data=d)
#' info.tempo <- vcov2(e.lm, bias.correct = TRUE)
#' info.tempo[names(coef(e.lm)),names(coef(e.lm))] - vcov(e.lm)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' vcov.tempo <- vcov2(e.lvm, bias.correct = FALSE)
#' round(vcov.tempo %*% information(e.lvm), 5)
#'
#' @concept small sample inference
#' @export
`information2` <-
  function(object, ...) UseMethod("information2")

## ** information2.lm
#' @rdname information2
#' @export
information2.lm <- function(object, ...){
    out <- .information2(dmu,
                         dOmega,
                         Omega,
                         n.corrected = object$sCorrect$n,
                         index.Omega = object$sCorrect$cluster$index.endogenous,
                         leverage,
                         n.cluster = object$sCorrect$cluster$n.cluster,,
                         grid.meanparam,
                         n.grid.meanparam,
                         grid.varparam,
                         n.grid.varparam,
                         name.param,
                         n.param)

    return(solve(vcov2(object, ...)))
}

## ** information2.gls
#' @rdname information2
#' @export
information2.gls <- information2.lm

## ** information2.lme
#' @rdname information2
#' @export
information2.lme <- information2.lm




## * vcov2
## ** vcov2 (documentation)
#' @title  Extract the Variance Covariance Matrix of the Model Parameters After Small Sample Correction
#' @description  Extract the variance covariance matrix of the model parameters from \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} objects.
#' Similar to \code{vcov} but with small sample correction (if any).
#' @name vcov2
#'
#' @param object a \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} object.
#' @param param [optional] the fitted parameters.
#' @param data [optional] the data set.
#' @param bias.correct [logical] should the standard errors of the coefficients be corrected for small sample bias? Only relevant if the \code{sCorrect} function has not yet be applied to the object.
#' @param ... arguments to be passed to \code{sCorrect}.
#' 
#' @details If argument \code{p} or \code{data} is not null, then the small sample size correction is recomputed to correct the influence function.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @return A matrix.
#' 
#' @examples
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#'
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- sequence.lvm(0)
#' set.seed(10)
#' d <- lava::sim(m,n)
#'
#' ## linear model
#' e.lm <- lm(formula.lvm,data=d)
#' vcov.tempo <- vcov2(e.lm, bias.correct = TRUE)
#' vcov.tempo[rownames(vcov(e.lm)),colnames(vcov(e.lm))]/vcov(e.lm)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' vcov.tempo <- vcov2(e.lvm, bias.correct = FALSE)
#' vcov.tempo/vcov(e.lvm)
#'
#' @concept small sample inference
#' @export
`vcov2` <-
  function(object, ...) UseMethod("vcov2")

## ** vcov2.lm
#' @rdname vcov2
#' @export
vcov2.lm <- function(object, param = NULL, data = NULL, ...){
    
    if(is.null(param) && is.null(data)){
        out <- vcov(object)
    }else{
        browser()
        ## sCorrect(object, param = param, data = data, df = FALSE, ...) <- bias.correct
    }
    
    return(out)
}

## ** vcov2.gls
#' @rdname vcov2
#' @export
vcov2.gls <- vcov2.lm

## ** vcov2.lme
#' @rdname vcov2
#' @export
vcov2.lme <- vcov2.lm

## ** vcov2.lvmfit
#' @rdname vcov2
#' @export
vcov2.lvmfit <- function(object, param = NULL, data = NULL, ...){

    if(is.null(param) && is.null(data)){
        out <- vcov(object)
    }else{
        out <- solve(information(object, p = param, data = data))
    }

    return(out)
}


## ** vcov2.sCorrect
#' @rdname vcov2
#' @export
vcov2.sCorrect <- function(object, param = NULL, data = NULL, ...){
    if(is.null(param) && is.null(data)){
        out <- object$sCorrect$vcov.param
    }else{
        out <- .info2vcov(object, p = param, data = data, attr.info = FALSE)
    }
}

## ** .info2vcov (helper)
#' @title Inverse the Information Matrix
#' @description Compute the inverse of the information matrix.
#' @name vcov2-internal
#'
#' @param attr.info [logical] should the information matrix be returned as an attribute?
#' @param ... arguments passed to .information2
#' 
#' @keywords internal
.info2vcov <- function(attr.info = FALSE, ...){
    iInfo <- .information2(...)
    iVcov <- try(chol2inv(chol(iInfo)), silent = TRUE)
    if(inherits(iVcov, "try-error")){
        iVcov <- solve(iInfo)
    }
    if(attr.info){
        attr(iVcov,"information") <- iInfo
    }
    return(iVcov)
}




######################################################################
### object-sCorrect.R ends here

### sCorrect-coef2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (10:14) 
## Version: 
## Last-Updated: feb 11 2020 (17:21) 
##           By: Brice Ozenne
##     Update #: 176
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Export Mean and Variance Coefficients
#' @description Export mean and variance coefficients from a \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit}.
#' Similar to \code{coef} but extract all model coefficients, including variance, correlation, and random effects,
#' possibly after application of the small sample correction.
#' @name coef2
#'
#' @param object a \code{lm}, \code{gls}, \code{lme}, \code{lvmfit} object.
#' @param as.lava [logical] Should the order and the name of the coefficient be the same as those obtained using coef with type = -1.
#' Only relevant for \code{lvmfit} objects.
#' 
#' @details For \code{gls} and \code{lme} models, the variance coefficients that are exported are the residual variance of each outcome. 
#' This is \eqn{\sigma^2} for the first one and \eqn{k^2 \sigma^2} for the remaining ones.
#'
#' @return A numeric vector named with the names of the coefficient with three attributes:
#' \itemize{
#' \item mean.coef: the name of the mean coefficients.
#' \item var.coef: the name of the variance coefficients.
#' \item cor.coef: the name of the correlation coefficients.
#' }
#'
#' @concept extractor
#' @keywords smallSampleCorrection
#' @export
`coef2` <-
    function(object, as.lava) UseMethod("coef2")


## * Examples
#' @rdname coef2
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
#' coef2(e.lm)
#'
#' #### gls model ####
#' e.gls1 <- gls(Y1~X1, data = dW)
#' coef2(e.gls1)
#' 
#' e.gls5 <- gls(Y~X1, weight = varIdent(form=~1|time2),
#'               correlation = corSymm(form=~time|id), data = dL)
#' coef2(e.gls5)
#'
#' #### lme model ####
#' e.lme1 <- lme(Y~X1, random = ~1|id, data = dL)
#' coef2(e.lme1)
#' 
#' e.lme2 <- lme(Y~X1, random = ~1|id, data = dL,
#'               weight = varIdent(form=~1|time2))
#' coef2(e.lme2)
#' 
#' #### latent variable models ####
#' e.lvm <- estimate(lvm(c(Y1,Y2,Y3) ~ 1*eta + X1, eta ~ Z1), data = dW)
#' coef2(e.lvm)
#' coef(e.lvm)

## * coef2.lm
#' @rdname coef2
#' @export
coef2.lm <- function(object, as.lava = TRUE){
    if(is.null(object$sCorrect)){
        out <- .coef2(object)
    }else{
        out <- object$sCorrect$param
    }
    return(out)
}

## * coef2.gls
#' @rdname coef2
#' @export
coef2.gls <- coef2.lm

## * coef2.lme
#' @rdname coef2
#' @export
coef2.lme <- coef2.lm

## * coef2.lvmfit
#' @rdname coef2
#' @export
coef2.lvmfit <- function(object, as.lava = TRUE){
    if(as.lava){
        out <- stats::coef(object)
    }else{
        out <- .coef2(object)
    }
    return(out)
}



## * coef2.sCorrect
#' @rdname coef2
#' @export
coef2.sCorrect <- function(object, as.lava = TRUE){
    out <- object$sCorrect$param
    if(as.lava == FALSE){ 
        out <- out[names(object$sCorrect$skeleton$originalLink2param)]
        names(out) <- as.character(object$sCorrect$skeleton$originalLink2param)
    }
    return(out)
}

## * .coef2
.coef2 <- function(object){

    if(inherits(object,"lm")){        
        coef.object <- coef(object)
        out <- c(coef.object,sigma2=mean(residuals(object)^2))
        attr(out, "mean.coef") <- names(coef.object)
        attr(out, "var.coef") <- "sigma2"

    }else if(inherits(object,"gls") || inherits(object,"lme")){

        ## ** mean coefficients
        if(inherits(object,"gls")){
            mean.coef <- stats::coef(object)
        }else if(inherits(object,"lme")){
            mean.coef <- nlme::fixef(object)
        }
            
        ## ** variance coefficients
        var.coef <- c(sigma2 = stats::sigma(object)^2)
        if(!is.null(object$modelStruct$varStruct)){
            var.coef <- c(var.coef,
                          stats::coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))          
        }

        ## ** random effect coefficients
        if(inherits(object,"lme")){
            random.coef <- as.double(nlme::getVarCov(object))    
            names(random.coef) <- paste0("ranCoef",1:length(random.coef))
        }else{
            random.coef <- NULL
        }
        
        ## ** correlation coefficients
        if(!is.null(object$modelStruct$corStruct)){
            cor.coef <- stats::coef(object$modelStruct$corStruct, unconstrained = FALSE)

            n.var <- length(var.coef)
            n.cor <- length(cor.coef)

        ## check unstructured covariance matrix
        if(!is.null(object$modelStruct$varStruct) && ((n.var*(n.var-1))/2 == n.cor)){

            vecgroup <- attr(unclass(object$modelStruct$corStruct), "group")
            veccov.cor <- unname(unlist(attr(object$modelStruct$corStruct, "covariate")))
            veccov.var <- attr(object$modelStruct$varStruct, "groups")

            table.covvar <- table(veccov.cor,veccov.var)
            newlevels.cor <- colnames(table.covvar)[apply(table.covvar, 1, which.max)]
            veccov.cor2 <- factor(veccov.cor, levels = 0:max(veccov.cor), labels = newlevels.cor)
            
            if(identical(as.character(veccov.cor2),as.character(veccov.var))){
                table.cor.coef <- utils::combn(newlevels.cor, m = 2)
                cor.coefName <- apply(table.cor.coef, MARGIN = 2, FUN = paste, collapse = "")
                names(cor.coef) <- paste0("corCoef",cor.coefName)

            }else{
                table.cor.coef <- NULL
                names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
            }
            
        }else{
            table.cor.coef <- NULL
            names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
        }
        
        
        }else{
            table.cor.coef <- NULL
            cor.coef <- NULL
        }
        out <- c(mean.coef, cor.coef, var.coef, random.coef)
        attr(out, "mean.coef") <- names(mean.coef)
        attr(out, "var.coef") <- names(var.coef)
        attr(out, "cor.coef") <- names(cor.coef)
        attr(out, "ran.coef") <- names(random.coef)
        attr(out, "table.cor.coef") <- table.cor.coef
        
    }else if(inherits(object,"lvmfit")){
        tempo <- stats::coef(object, type = 2, labels = 1)
        out <- setNames(tempo[,1],rownames(tempo))
        attr(out, "mean.coef") <- rownames(tempo)[attr(tempo,"type")!="variance"]
        attr(out, "var.coef") <- rownames(tempo)[attr(tempo,"type")=="variance"]
    }

    return(out)
}

######################################################################
### sCorrect-coef2.R ends here

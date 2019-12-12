### sCorrect-coef2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (10:14) 
## Version: 
## Last-Updated: dec 11 2019 (16:22) 
##           By: Brice Ozenne
##     Update #: 44
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
`coef2` <-
    function(object, ssc, labels) UseMethod("coef2")


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

## * coef2.lm
#' @rdname coef2
coef2.lm <- function(object, ssc = TRUE, labels){
    if(ssc){
        if(is.null(object$sCorrect) || (object$sCorrect$ssc != ssc)){
            object <- sCorrect(object, ssc = ssc, df = FALSE)
        }
        out <- object$sCorrect$coef
    }else{
        coef.object <- coef(object)
        out <- c(coef.object,sigma2=sigma(object)^2)
        attr(out, "mean.coef") <- names(coef.object)
        attr(out, "var.coef") <- "sigma2"
    }

    return(out)
}

## * coef2.gls
#' @rdname coef2
coef2.gls <- function(object, ssc = TRUE, labels){
    if(ssc){
        if(is.null(object$sCorrect) || (object$sCorrect$ssc != ssc)){
            object <- sCorrect(object, ssc = ssc, df = FALSE)
        }
        out <- object$sCorrect$coef
    }else{
        
    ## *** mean coefficients
    mean.coef <- stats::coef(object)

    ## *** variance coefficients
    var.coef <- c(sigma2 = stats::sigma(object)^2)
    if(!is.null(object$modelStruct$varStruct)){
        var.coef <- c(var.coef,
                      stats::coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))          
    }

    ## *** correlation coefficients
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

                cor.coefName <- apply(utils::combn(newlevels.cor, m = 2), MARGIN = 2, FUN = paste, collapse = "")
                names(cor.coef) <- paste0("corCoef",cor.coefName)

            }else{
                names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
            }
            
            
        }else{
            names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
        }
        
        
    }else{
        cor.coef <- NULL
    }

        out <- c(mean.coef, cor.coef, var.coef)
        attr(out, "mean.coef") <- names(mean.coef)
        attr(out, "var.coef") <- names(var.coef)
        attr(out, "cor.coef") <- names(cor.coef)
    }
    return(out)
}




## * coef2.lme
#' @rdname coef2
coef2.lme <- function(object, ssc = TRUE, labels){

    if(ssc){
        if(is.null(object$sCorrect) || (object$sCorrect$ssc != ssc)){
            object <- sCorrect(object, ssc = ssc, df = FALSE)
        }
        out <- object$sCorrect$coef
    }else{
        ## *** mean coefficients
        mean.coef <- nlme::fixef(object)

    ## *** variance coefficients
    var.coef <- c(sigma2 = stats::sigma(object)^2)
    if(!is.null(object$modelStruct$varStruct)){
        var.coef <- c(var.coef,
                      stats::coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))   
    }

    ## *** random effect coefficients
    random.coef <- as.double(nlme::getVarCov(object))    
    names(random.coef) <- paste0("ranCoef",1:length(random.coef))

    ## *** correlation coefficients
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

                cor.coefName <- apply(utils::combn(newlevels.cor, m = 2), MARGIN = 2, FUN = paste, collapse = "")
                names(cor.coef) <- paste0("corCoef",cor.coefName)

            }else{
                names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
            }
            
            
        }else{
            names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
        }
        
        
        }else{
            cor.coef <- NULL
        }
    
    out <- c(mean.coef, cor.coef, var.coef, random.coef)

            attr(out, "mean.coef") <- names(mean.coef)
            attr(out, "var.coef") <- names(var.coef)
            attr(out, "cor.coef") <- names(cor.coef)
            attr(out, "ran.coef") <- names(random.coef)
            return(out)
        }
}

## * coef2.lvmfit
#' @rdname coef2
coef2.lvmfit <- function(object, ssc = TRUE, labels = lava.options()$coef.names){
    if(ssc){
        if(is.null(object$sCorrect) || (object$sCorrect$ssc != ssc)){
            object <- sCorrect(object, ssc = ssc, df = FALSE)
        }
        out <- object$sCorrect$coef
    }else{
        tempo <- coef(object, type = 2, labels = labels)
        out <- tempo[,1]
        attr(out, "mean.coef") <- rownames(tempo)[attr(tempo,"type")!="variance"]
        attr(out, "var.coef") <- rownames(tempo)[attr(tempo,"type")=="variance"]
        attr(out, "cor.coef") <- NULL
        return(out)
    }
}

## * coef.sCorrect
#' @rdname coef2
coef.sCorrect <- function(object, ssc = TRUE, labels = lava.options()$coef.names){
    class(object) <- setdiff(class(object),"sCorrect")
    return(coef2(object, ssc = ssc, labels = labels))

}

######################################################################
### sCorrect-coef2.R ends here

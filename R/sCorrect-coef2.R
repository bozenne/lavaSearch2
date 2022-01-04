### sCorrect-coef2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (10:14) 
## Version: 
## Last-Updated: Jan  4 2022 (17:45) 
##           By: Brice Ozenne
##     Update #: 261
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Model Coefficients With Small Sample Correction
#' @description Extract the coefficients from the latent variable model.
#' Similar to \code{lava::compare} but with small sample correction.
#' @name coef2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param as.lava [logical] if \code{TRUE}, uses the same names as when using \code{stats::coef}.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients (\code{"none"}, \code{"residual"}, \code{"cox"}). Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object. 
#'
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the model coefficients.
#' 
#' @return A numeric vector named with the names of the coefficients.
#' 
#' @seealso \code{\link{estimate2}} to obtain \code{lvmfit2} objects.
#' 
#' @concept extractor
#' @keywords smallSampleCorrection
#' @export
`coef2` <-
    function(object, as.lava, ...) UseMethod("coef2")


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
#' #### latent variable models ####
#' e.lvm <- estimate(lvm(c(Y1,Y2,Y3) ~ 1*eta + X1, eta ~ Z1), data = dW)
#' coef(e.lvm)
#' coef2(e.lvm)
#' coef2(e.lvm, as.lava = FALSE)

## * coef2.lvmfit
#' @rdname coef2
coef2.lvmfit <- function(object, as.lava = TRUE, ssc = lava.options()$ssc, ...){

    return(coef(estimate2(object, ssc = ssc, ...), as.lava = as.lava))

}

## * coef2.lvmfit2
#' @rdname coef2
coef2.lvmfit2 <- function(object, as.lava = TRUE, ...){

    dots <- list(...)
    if(length(dots)>1){ ## for the print function

        ## extract structure from lava
        object0 <- object
        class(object0) <- setdiff(class(object), "lvmfit2")
        out <- do.call(stats::coef, args = c(list(object0),dots))
        ## get the full name from lava
        dots$symbol <- NULL
        out.names <- rownames(do.call(stats::coef, args = c(list(object0),dots)))
        ## extract value from lavaSearch2
        res <- confint(object, as.lava = FALSE)[names(object$sCorrect$skeleton$originalLink2param),]
        rownames(res) <- as.character(object$sCorrect$skeleton$originalLink2param)

        out[,"Estimate"] <- res[match(out.names,object$sCorrect$skeleton$originalLink2param),"estimate"]
        out[,"Std. Error"] <- res[match(out.names,object$sCorrect$skeleton$originalLink2param),"std.error"]
        out[,"Z-value"] <- res[match(out.names,object$sCorrect$skeleton$originalLink2param),"statistic"]
        colnames(out)[colnames(out)=="Z-value"] <- "t-value"
        out[,"P-value"] <- res[match(out.names,object$sCorrect$skeleton$originalLink2param),"p.value"]
    }else{        
        out <- object$sCorrect$param
        if(as.lava == FALSE){
            out <- out[names(object$sCorrect$skeleton$originalLink2param)]
            names(out) <- as.character(object$sCorrect$skeleton$originalLink2param)
        }
    }
    return(out)
}

## * coef.lvmfit2
#' @rdname coef2
#' @export
coef.lvmfit2 <- coef2.lvmfit2

######################################################################
### sCorrect-coef2.R ends here

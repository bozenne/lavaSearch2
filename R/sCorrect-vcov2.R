### sCorrect-vcov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 11 2019 (13:55) 
## Version: 
## Last-Updated: dec 13 2019 (10:44) 
##           By: Brice Ozenne
##     Update #: 8
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

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
  function(object, param, data, ssc) UseMethod("vcov2")

## * vcov2.lm
#' @rdname vcov2
#' @export
vcov2.lm <- function(object, param = NULL, data = NULL, ssc = TRUE){
    
    if(is.null(object$sCorrect) || !is.null(param) || !is.null(data) || (object$sCorrect$ssc != ssc)){
        object <- sCorrect(object, param = param, data = data, ssc = ssc, df = FALSE)
    }
    
    return(.info2vcov(object$sCorrect$information))
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
vcov2.lvmfit <- vcov2.lm


## ** vcov2.sCorrect
#' @rdname vcov2
#' @export
vcov2.sCorrect <- function(object, param = param, data = data, ssc = ssc){
    class(object) <- setdiff(class(object),"sCorrect")
    return(information2(object, param = param, data = data, ssc = ssc))
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
.info2vcov <- function(information, attr.info = FALSE){
    vcov <- try(chol2inv(chol(information)), silent = TRUE)
    if(inherits(vcov, "try-error")){
        vcov <- try(solve(information), silent = FALSE)
    }
    if(attr.info){
        attr(vcov,"information") <- information
    }
    if(!inherits(vcov, "try-error")){
        dimnames(vcov) <- dimnames(information)
    }
    return(vcov)
}


######################################################################
### sCorrect-vcov2.R ends here

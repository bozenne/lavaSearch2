### method-sCorrect.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  2 2019 (13:39) 
## Version: 
## Last-Updated: nov 18 2019 (10:09) 
##           By: Brice Ozenne
##     Update #: 6
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * information2.lm2
#' @rdname information2
#' @export
information2.lm2 <- function(object, ...){
    return(object)
}

## * information2.gls2
#' @rdname information2
#' @export
information2.gls2 <- information2.lm2

## * information2.lme2
#' @rdname information2
#' @export
information2.lme2 <- information2.lm2

## * information2.lvmfit
#' @rdname information2
#' @export
information2.lvmfit2 <- information2.lm2

## * compare2 
## * effects2 
## * summary2 

## * nobs
nobs.gls2 <- function(object){
    return(NROW(object$sCorrect$score))
}
nobs.lme2 <- nobs.gls2


######################################################################
### method-sCorrect.R ends here

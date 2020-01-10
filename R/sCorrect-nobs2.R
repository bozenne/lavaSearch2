### method-sCorrect.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  2 2019 (13:39) 
## Version: 
## Last-Updated: jan  8 2020 (11:11) 
##           By: Brice Ozenne
##     Update #: 17
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * nobs2
#' @export
`nobs2` <-
    function(object) UseMethod("nobs2")

## * nobs2.lm
#' @rdname nobs2
#' @export
nobs2.lm <- function(object){
    if(!is.null(object$sCorrect)){
        out <- object$sCorrect$cluster$n.cluster
    }else{
        data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                            envir = parent.env(environment()), rm.na = TRUE)
        endogenous <- endogenous(object)

        out <- .getGroups2(object, data = data, endogenous = endogenous)$n.cluster
    }
    return(out)
}

## * nobs2.gls
#' @rdname nobs2
#' @export
nobs2.gls <- nobs2.lm

## * nobs2.lme
#' @rdname nobs2
#' @export
nobs2.lme <- nobs2.lm

## * nobs2.lvmfit
#' @rdname nobs2
#' @export
nobs2.lvmfit <- nobs2.lm

######################################################################
### method-sCorrect.R ends here

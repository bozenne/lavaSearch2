### df.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 23 2017 (12:27) 
## Version: 
## last-updated: jun 27 2017 (11:38) 
##           By: Brice Ozenne
##     Update #: 21
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Approximate degree of freedoms of a model
#' @description Approximate the degree of freedoms of a model
#' @name df.residual
#' 
#' @param object model object
#' @param conservative If true the total number of parameter is substracted from the number of observations. If false only the number of mean parameters is substracted.
#' @param ... additional arguments
#' @examples
#' n <- 20
#' 
#' set.seed(10)
#' m <- lvm(Y~X1+X2+X3)
#' d <- sim(m, n)
#' e <- estimate(m, data = d)
#' df.residual(e, conservative = FALSE)
#' df.residual(e, conservative = TRUE)

#' @rdname df.residual
#' @method df.residual lvmfit
#' @export
df.residual.lvmfit <- function(object, conservative,...) {

    n <- object$data$n
    object.coef <- coef(object, level = 9)
    p.mean <- sum(attr(object.coef,"type")=="regression")
    p.intercept <- sum(attr(object.coef,"type")=="intercept")
    p.variance <- sum(attr(object.coef,"type")=="variance")

    if(conservative){
        p.effective <- p.mean+p.intercept
    }else{
        p.effective <- p.mean+p.intercept+p.variance
    }

    
    return(n-p.effective)
}

#' @rdname df.residual
#' @method df.residual coxph
#' @export
df.residual.coxph <- function(object, ...) {
    n <- riskRegression::CoxN(object)    
    object.coef <- coef(object)
    p.effective <- length(object.coef)
    
    return(n-p.effective)
}
#----------------------------------------------------------------------
### df.R ends here

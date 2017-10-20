### df.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 23 2017 (12:27) 
## Version: 
## last-updated: okt  5 2017 (11:32) 
##           By: Brice Ozenne
##     Update #: 25
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
#' m <- lvm(Y~X1+X2+X3+X4)
#' d <- lava::sim(m, n)
#' e <- estimate(m, data = d)
#' df.residual(e, conservative = FALSE)
#' df.residual(e, conservative = TRUE)
#'
#' 
#' @rdname df.residual
#' @method df.residual lvmfit
#' @export
df.residual.lvmfit <- function(object, conservative,...) {

    n <- object$data$n
    object.coef <- coef(object, level = 9)
    index.keep <- which(rownames(object.coef) %in% names(coef(object)))
    p.mean <- sum(attr(object.coef,"type")[index.keep]=="regression")
    p.intercept <- sum(attr(object.coef,"type")[index.keep]=="intercept")
    p.variance <- sum(attr(object.coef,"type")[index.keep]=="variance")

    if(conservative){
      p.effective <- p.mean+p.intercept+p.variance
    }else{
      p.effective <- p.mean+p.intercept
    }

    
    return(n-p.effective)
}

#' @rdname df.residual
#' @method df.residual coxph
#' @export
df.residual.coxph <- function(object, ...) {
    n <- riskRegression::coxN(object)    
    object.coef <- coef(object)
    p.effective <- length(object.coef)
    
    return(n-p.effective)
}
#----------------------------------------------------------------------
### df.R ends here

### df.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 23 2017 (12:27) 
## Version: 
## last-updated: jun 26 2017 (10:36) 
##           By: Brice Ozenne
##     Update #: 12
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Approximate the degree of freedoms of a lvm model
#' @description Approximate the degree of freedoms of a lvm model
#'
#' @name degreeOfFreedom
#' 
#' @param x model object
#' @param conservative If true the total number of parameter is substracted from the number of observations. If false only the number of mean parameters is substracted.
#' @param ... additional arguments
#' @examples
#' n <- 20
#' 
#' set.seed(10)
#' m <- lvm(Y~X1+X2+X3)
#' d <- sim(m, n)
#' e <- estimate(m, data = d)
#' degreeOfFreedom(e, conservative = FALSE)
#' degreeOfFreedom(e, conservative = TRUE)
#' 
#' @export
degreeOfFreedom <- function(x,...) UseMethod("degreeOfFreedom")

#' @rdname degreeOfFreedom
#' @export
degreeOfFreedom.lvmfit <- function(x,conservative,...) {

    n <- x$data$n
    x.coef <- coef(x, level = 9)
    p.mean <- sum(attr(x.coef,"type")=="regression")
    p.intercept <- sum(attr(x.coef,"type")=="intercept")
    p.variance <- sum(attr(x.coef,"type")=="variance")

    if(conservative){
        p.effective <- p.mean+p.intercept
    }else{
        p.effective <- p.mean+p.intercept+p.variance
    }

    
    return(n-p.effective)
}
#----------------------------------------------------------------------
### df.R ends here

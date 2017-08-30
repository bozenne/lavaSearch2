### qqplot2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 30 2017 (09:26) 
## Version: 
## last-updated: aug 30 2017 (11:08) 
##           By: Brice Ozenne
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ documentation

#' @title QQplot for the residuals of a lvm object
#' @description QQplot for the residuals of a lvm object
#' 
#' @name qqplot2
#'
#' @param object a lvm model.
#' @param variables the variable for which the residuals should be displayed.
#' @param mfrow how to divide the window. See \code{mfrow}.
#' @param type the function used to display the qqplot. Can be qqtest or qqnorm.
#' @param centralPercents argument passed to \code{qqtest}. See the help of \code{\link{qqtest}}.
#' @param ... additional arguments to be passed to qqtest.
#' 
#' @details 
#' Simulation is based on a multivariate truncated normal law (even though it is not satifying for the variance components)
#' 
#' @return a data frame/cvlvm object containing the convergence status (by default 0 indicates successful convergence, see ?optim), the value of the log-likelihood and the estimated parameters (in columns) for each initialization (in rows)
#' 
#' @examples 
#' m <- lvm(list(y~v1+v2+v3+v4,c(v1,v2,v3,v4)~x))
#' covariance(m) <- v1~v2+v3+v4
#' latent(m) <- ~ x
#' dd <- sim(m,100) ## Simulate 100 observations from model
#' e <- estimate(m, dd) ## Estimate parameters
#'
#' qqplot2(e)
#' 
#' @export
qqplot2 <- function (object, ...) {
  UseMethod("qqplot2", object)
}
# }}}

#' @rdname qqplot2
#' @export
qqplot2.lvmfit <- function(object, variables = NULL, mfrow = NULL,
                             type = "qqtest",  centralPercents = 0.95,...){

    M.res <- predict(object, residual = TRUE)
    name.vars <- colnames(M.res)
    
    if(!is.null(variables)){
        if(any(variables %in% name.vars == FALSE)){
            stop("unknown variable(s): ",paste(variables[variables %in% name.vars == FALSE], collapse = " "),"\n",
                 "endogenous variables: ",paste(endogenous(object), collapse = " "),"\n",
                 "latent variables: ",paste(latent(object), collapse = " "),"\n")
        }
        M.res <- M.res[,variables,drop=FALSE]
        name.vars <- variables
    }
    if(type %in% c("qqtest","qqnorm") == FALSE){
        stop("wrong specification of type \n",
             "must be \"qqtest\" or \"qqnorm\" \n")
    }

    n.var <- NCOL(M.res)       
    if(is.null(mfrow)){
        mfrow <- c(round(sqrt(n.var)), ceiling(n.var/round(sqrt(n.var))))
    }
    op <- par(mfrow = mfrow)
    sapply(1:n.var, function(row){
        resid <- na.omit(M.res[,row])
        main <- name.vars[row]
        if(all(resid < 1e-5)){
            plot(0,0, col = "white", axes = FALSE, xlab = "", ylab = "", main = main)
            text(0,0,"all residuals < 1e-5")
        }else if(type == "qqtest"){
            qqtest::qqtest(resid, main = name.vars[row],
                           centralPercents = centralPercents,
                           ...)
        }else if(type == "qqnorm"){
            qqnorm(M.res[,row], main = name.vars[row])
        }
    })
    par(op)

    return(invisible(M.res))
}

#----------------------------------------------------------------------
### qqplot2.R ends here

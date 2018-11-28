### clean.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 27 2018 (14:35) 
## Version: 
## Last-Updated: nov 28 2018 (11:32) 
##           By: Brice Ozenne
##     Update #: 7
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * clean (Documentation)
#' @title Simplify a lvm object
#' @description Remove variables with no link.
#' @name clean
#' 
#' @param x \code{lvm}-object
#' @param rm.exo  should the exogenous variables with no links be removed from the object ? 
#' @param rm.endo  should the endogenous variables with no links be removed from the object ?
#' @param rm.latent  should the latent variables with no links be removed from the object ?
#' @param ... additional arguments to lower level functions
#'
#' @export
`clean` <-
  function(x,...) UseMethod("clean")

## * clean (examples)
#' @rdname clean
#' @examples 
#' 
#' m <- lvm()
#' m <- regression(m, x=paste0("x",1:5),y="y1")
#' m <- regression(m, x=paste0("x",1:5),y="y2")
#' covariance(m) <- y1~y2
#'
#' cancel(m) <- y1 ~ x1
#' cancel(m) <- y2 ~ x1
#' clean(m)
#'
#' m <- lvm(y1 ~ eta + x1, y2 ~ eta, y3 ~ eta + x2)
#' latent(m) <- ~eta
#' clean(m)
#' m
#' cancel(m) <- y1 ~ eta
#' cancel(m) <- y2 ~ eta
#' cancel(m) <- y3 ~ eta
#' clean(m)

## * clean.lvm
#' @rdname clean
#' @export
clean.lvm <- function(x, rm.exo = TRUE, rm.endo = TRUE, rm.latent = TRUE, ...){

    myhooks <- gethook_lavaSearch2("clean.hooks")
    for (f in myhooks) {
        res <- do.call(f, list(x=x, rm.exo=rm.exo,...))
        if("x" %in% names(res)){ x <- res$x }
        if("rm.exo" %in% names(res)){ rm.exo <- res$rm.exo }
        if("rm.endo" %in% names(res)){ rm.endo <- res$rm.endo }
    }  
    
    var.latent <- latent(x)
    var.exogenous <- exogenous(x)
    var.endogenous <- endogenous(x)

    M.reg <- x$index$A
    M.cov <- x$index$P - diag(diag(x$index$P))
    
    varClean <- NULL
    if(rm.exo && length(var.exogenous) > 0){
        indexClean.reg <- which(rowSums(M.reg[var.exogenous,,drop = FALSE]!=0)==0)
        indexClean.cov <- which(rowSums(M.cov[var.exogenous,,drop = FALSE]!=0)==0)
        indexClean <- intersect(indexClean.reg, indexClean.cov)
        varClean <- c(varClean, var.exogenous[indexClean])
    }
    if(rm.endo && length(var.endogenous)>0){
        indexClean.reg <- which(colSums(M.reg[,var.endogenous,drop = FALSE]!=0)==0)
        indexClean.cov <- which(colSums(M.cov[,var.endogenous,drop = FALSE]!=0)==0)
        indexClean <- intersect(indexClean.reg, indexClean.cov)
        varClean <- c(varClean, var.endogenous[indexClean])
    }
    if(rm.latent && length(var.latent)>0){
        indexClean.Rreg <- which(rowSums(M.reg[var.latent,,drop = FALSE]!=0)==0)
        indexClean.Rcov <- which(rowSums(M.cov[var.latent,,drop = FALSE]!=0)==0)
        indexClean.Creg <- which(colSums(M.reg[,var.latent,drop = FALSE]!=0)==0)
        indexClean.Ccov <- which(colSums(M.cov[,var.latent,drop = FALSE]!=0)==0)
        indexClean <- intersect(intersect(indexClean.Rreg, indexClean.Rcov),
                                intersect(indexClean.Creg, indexClean.Ccov))
        varClean <- c(varClean, var.latent[indexClean])
    }
    
    if(length(varClean)>0){
        x <- kill(x, value =  varClean, ...)
    }
    return(x)
}


######################################################################
### clean.R ends here

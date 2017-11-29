### mlf2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2017 (12:56) 
## Version: 
## Last-Updated: nov 29 2017 (17:31) 
##           By: Brice Ozenne
##     Update #: 32
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * mlf2
#' @title Simultaneous Inference for Multiple Models
#' @description Calculation of correlation between test statistics from multiple marginal models using the score decomposition. Similar to \code{multcomp::mmm} but with small sample corrections.
#'
#' @param .. A names argument list containing the fitted model.
#' See the documentation of \code{multcomp::mmm} for more details.
#'
#' @export
mlf2 <- function(...) {

     ret <- list(...)
     class(ret) <- "mlf2"
     ret
}

## * glht.mlf2
## copy of multcomp:::glht.mlf, allowing mmm2 objects
glht.mlf2 <- function(model, linfct, ...) {

    stopifnot(inherits(model, c("mmm2")))
    if (length(linfct) == 1) {
        linfct <- linfct[rep(1, length(model))]
        names(linfct) <- names(model)
    }
    if (!isTRUE(all.equal(names(model), names(linfct))))
        stop("names of ", sQuote("model"), " and ", sQuote("linfct"),
             " are not identical")

    K <- lapply(names(model), function(i) 
        glht(model[[i]], linfct = linfct[[i]], ...)$linfct)
    for (k in 1:length(K)) {
        rownames(K[[k]]) <- paste(names(model)[k], 
                                  rownames(K[[k]]), sep = ": ")
        colnames(K[[k]]) <- paste(names(model)[k], 
                                  colnames(K[[k]]), sep = ": ")
    }
    ## do not add variance in the first run because cannot pass addtional arguments
    out <- multcomp_glht.matrix(model, linfct = multcomp_.bdiag(K), ...)
    
    ## add variance
    out$vcov <- vcov(model, return.null = FALSE, ...)

    ## add degrees of freedom
    df <- unlist(lapply(model,dfVariance,...))
    out$df <- as.double(round(median(df)))
    return(out)
}

##----------------------------------------------------------------------
### mlf2.R ends here

### mlf2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2017 (12:56) 
## Version: 
## Last-Updated: jan  9 2018 (17:57) 
##           By: Brice Ozenne
##     Update #: 47
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
#' @param ... A names argument list containing the fitted model.
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

    ## ** check
    stopifnot(inherits(model, c("mmm2")))
    if (length(linfct) == 1) {
        linfct <- linfct[rep(1, length(model))]
        names(linfct) <- names(model)
    }
    if (!isTRUE(all.equal(names(model), names(linfct))))
        stop("names of ", sQuote("model"), " and ", sQuote("linfct"),
             " are not identical")

    ## ** define contrast matrix
    n.model <- length(model)
    name.model <- names(model)
    if(is.null(name.model)){
        name.model <- 1:n.coef
        names(model) <- name.model
    }
    
    K <- lapply(names(model), function(i) {
        glht(model[[i]], linfct = linfct[[i]], ...)$linfct
    })
    
    for (iK in 1:n.model) {
        rownames(K[[iK]]) <- paste(names(model)[iK], rownames(K[[iK]]), sep = ": ")
        colnames(K[[iK]]) <- paste(names(model)[iK], colnames(K[[iK]]), sep = ": ")
    }
    
    ## ** do not add variance in the first run because cannot pass addtional arguments
    out <- multcomp_glht.matrix(model, linfct = multcomp_.bdiag(K), ...)
    
    ## ** add variance
    ## call vcov.mmm2
    out$vcov <- vcov(model, return.null = FALSE, ...)

    ## add degrees of freedom
    names(K) <- names(model)
    df <- sapply(names(model), function(iName){ ## iName <- names(model)[1]
        if(identical(class(model[[iName]]),"lm") && "sigma2" %in% names(K[[iName]]) == FALSE){
            K[[iName]] <- cbind(K[[iName]], sigma2 = 0)
        }
        lTest(model[[iName]], C = K[[iName]], Ftest = FALSE, ...)$df
    })
    out$df <- as.double(round(median(df)))
    return(out)
}

##----------------------------------------------------------------------
### mlf2.R ends here

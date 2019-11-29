### formula.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 25 2019 (09:39) 
## Version: 
## Last-Updated: nov 25 2019 (10:07) 
##           By: Brice Ozenne
##     Update #: 13
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

formula.varStruct <- function(x, ...){
    return(attr(x, "formula"))
}
formula.corStruct <- function(x, ...){
    return(attr(x, "formula"))
}
formula.reStruct <- function(x, ...){
    n.random <- length(x)
    group.random <- names(x)
    ls.formula <- lapply(1:n.random, function(iN){ ## iN <- names(x)[1]
        as.formula(paste0(deparse(attr(x[[iN]],"formula")),"|", group.random[iN]))
    })
    return(ls.formula)
}


######################################################################
### formula.R ends here

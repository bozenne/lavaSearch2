### calcDDF.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (09:29) 
## Version: 
## last-updated: nov  9 2017 (18:19) 
##           By: Brice Ozenne
##     Update #: 53
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - calcDDF
#' @title  Compute the degree of freedom of the variance parameters
#' @description Compute the degree of freedom of the variance parameters
#' @name calcDDF
#' @export
`calcDDF` <-
  function(x, ...) UseMethod("calcDDF")

## * calcDDF.gls
#' @rdname calcDDF
#' @export
calcDDF.gls <- function(x, p, iid0,
                        adjust.residuals, power, ...){
    browser()
    x$modelStruct$reStruct
    
    n.coef <- NCOL(iid0)
    calcSigma <- function(iP){ # iP <- p
        Sigma.tempo <- crossprod(iid2(x = x, p = iP,
                                      adjust.residuals = adjust.residuals, power = power,
                                      return.df = FALSE, check.score = FALSE, ...))
        return(Sigma.tempo)             
    }
    args(iid2.lme)
    browser()
    
    calcSigma(p)
    Sigma.beta <- crossprod(iid0)

    numerator <- 2*Sigma.beta^2
    jacob.iid0 <- numDeriv::jacobian(func = calcSigma, x = p)
    browser()
    denominator <- matrix(rowSums((jacob.iid0 %*% Sigma.beta) * jacob.iid0), nrow = n.coef, ncol = n.coef)

    numerator / denominator
    ## jacob.iid0 %*% Sigma.beta %*% t(jacob.iid0)

        ## calcSigma <- function(iP){ # iP <- p
        ##     Sigma.tempo <- crossprod(iid2.lvmfit(x, p = iP, data = data, adjust.residuals = adjust.residuals, power = power,
        ##                                          use.information = use.information, Dmethod = Dmethod, return.df = FALSE,
        ##                                          check.score = FALSE, ...))
        ##     return(Sigma.tempo[2,2])             
        ## }
 
}

## * calcDDF.lme
#' @rdname calcDDF
#' @export
calcDDF.lme <- function(x, p, iid0,
                        adjust.residuals, power, ...){

### ** normalize arguments
    if(!is.null(x$modelStruct$corStruct)){
        stop("cannot handle lme objects when corStruct is not null \n")
    }
    if(length(getVarCov(x))>1){
        stop("cannot handle lme objects with more than one random effect \n")
    }
    if(!is.null(x$modelStruct$varStruct) && "varIdent" %in% class(x$modelStruct$varStruct) == FALSE){
        stop("can only handle varIdent class for the variance structure \n")
    }
### ** prepare
    name.fixef <- names(fixef(x))
    test.weigth <- !is.null(x$modelStruct$varStruct)
    if(test.weigth){
        name.weight <- paste0("weight",1:length(x$modelStruct$varStruct))
        vec.weight <- setNames(as.double(x$modelStruct$varStruct),name.weight)
    }else{
        vec.weight <- NULL
    }
    vec.allCoef <- c(fixef(x),
                     sigma2=x$sigma,
                     tau=as.double(getVarCov(x)),
                     vec.weight)

    calcSigma <- function(iP){ # iP <- vec.allCoef
        ## print(iP-vec.allCoef)
        
        ## update variance parameters in the model
        x$sigma <- as.double(iP["sigma2"])
        x$modelStruct$reStruct[[1]][] <- iP["tau"]/x$sigma^2
        x$modelStruct$varStruct[] <- iP[name.weight]        
        old2new <- 1/coef(x$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)
        attr(x$modelStruct$varStruct,"weights") <- old2new[attr(x$modelStruct$varStruct, "groups")]
           
        Sigma.tempo <- crossprod(iid2(x = x, p = iP[name.fixef],
                                      adjust.residuals = adjust.residuals, power = power,
                                      return.df = FALSE, check.score = FALSE, ...))

        ## print(Sigma.tempo)
        return(as.double(Sigma.tempo))             
    }
                                       
    Sigma.beta <- crossprod(iid0)
    ## Sigma.beta-calcSigma(vec.allCoef)
    ## vec.allCoef2 <- vec.allCoef
    ## vec.allCoef2[1] <- vec.allCoef2[1] + 1
    ## calcSigma(vec.allCoef2)-calcSigma(vec.allCoef)
    numerator <- 2*Sigma.beta^2
    browser()
    jacob.iid0 <- numDeriv::jacobian(func = calcSigma, x = vec.allCoef)
    denominator <- matrix(rowSums((jacob.iid0 %*% Sigma.beta) * jacob.iid0), nrow = n.coef, ncol = n.coef)

    numerator / denominator
    ## jacob.iid0 %*% Sigma.beta %*% t(jacob.iid0)

        ## calcSigma <- function(iP){ # iP <- p
        ##     Sigma.tempo <- crossprod(iid2.lvmfit(x, p = iP, data = data, adjust.residuals = adjust.residuals, power = power,
        ##                                          use.information = use.information, Dmethod = Dmethod, return.df = FALSE,
        ##                                          check.score = FALSE, ...))
        ##     return(Sigma.tempo[2,2])             
        ## }
 
}

#----------------------------------------------------------------------
### calcDDF.R ends here

## * calcDDF.lvmfit
#' @rdname calcDDF
#' @export
calcDDF.lvmfit <- function(x, p, data, iid0,
                           adjust.residuals, power, as.clubSandwich, ...){

### ** prepare
    
    calcSigma <- function(iP){ # iP <- vec.allCoef
        
        Sigma.tempo <- crossprod(iid2(x = x, p = iP, data = data,
                                      adjust.residuals = adjust.residuals, power = power, as.clubSandwich = as.clubSandwich,
                                      return.df = FALSE, check.score = FALSE, ...))

        ## print(Sigma.tempo)
        ## return(as.double(Sigma.tempo))
        return(diag(Sigma.tempo))             
    }

    calcSigma(p+c(0,0,0,0,1))
    Sigma.beta <- crossprod(iid0)
    numerator <- 2*Sigma.beta^2
    browser()
    jacob.iid0 <- numDeriv::jacobian(func = calcSigma, x = p)
    denominator <- matrix(rowSums((jacob.iid0 %*% Sigma.beta) * jacob.iid0), nrow = n.coef, ncol = n.coef)

    numerator / denominator
    ## jacob.iid0 %*% Sigma.beta %*% t(jacob.iid0)

        ## calcSigma <- function(iP){ # iP <- p
        ##     Sigma.tempo <- crossprod(iid2.lvmfit(x, p = iP, data = data, adjust.residuals = adjust.residuals, power = power,
        ##                                          use.information = use.information, Dmethod = Dmethod, return.df = FALSE,
        ##                                          check.score = FALSE, ...))
        ##     return(Sigma.tempo[2,2])             
        ## }
 
}

#----------------------------------------------------------------------
### calcDDF.R ends here

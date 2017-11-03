### iid2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (13:16) 
## Version: 
## last-updated: okt 27 2017 (10:00) 
##           By: Brice Ozenne
##     Update #: 306
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - iid2
#' @title  Extract i.i.d. decomposition from linear and latent variable models
#' @description  Extract i.i.d. decomposition from linear and latent variable models using leverage adjusted fitted residuals.
#' @name iid2
#'
#' @param x a linear model or a latent variable model
#' @param data [optional] data set.
#' @param indiv Should the score relative to each observation be exported? Otherwise the total score (i.e. sum over all observations) will be exported.#' 
#' @param adjust.residuals Should the leverage-adjusted residuals be used to compute the score? Otherwise the raw residuals will be used.#' 
#' @param power the exponent used for computing the leverage-adjusted residuals. See the details section.
#' @param return.df Should the degree of freedom be computed?
#' @param ... arguments to be passed to \code{score2}.
#'
#' @details
#' Leverage-adjusted residuals have been shown to improve the coverage of robust standard errors in small samples.
#' They are computed according to the formula:
#' \eqn{e_adj = \frac{e}{(1-h_{ii})^\alpha}}
#'
#' \code{power = 0} (i.e. code{adjust.residuals = FALSE}) corresponds to the use of the raw residuals.
#' It can be shown that the variance-covariance estimator is biased downward for the mean parameters in LVM.
#' 
#' \code{power=0.5} corrects for the bias of the variance-covariance estimator,
#' at least for the mean parameters (Kauermann and Carroll, 2001)..
#' This corresponds to the correction "CR2" in the clubSandwich package.
#' 
#' \code{power=1} approximates the variance-covariance jackknife estimator (Bell and McCaffrey, 2002).
#' This corresponds to the correction "CR3" in the clubSandwich package.
#' 
#' @references
#' Bell, R. M., & McCaffrey, D. F. Bias reduction in standard errors for linear regression with multi-stage samples. Survey Methodology, 28(2), 169-181 (2002). \cr
#' Kauermann G. and Carroll R. J. A note on the efficiency of sandwich covariance matrix estimation. Journal of the American Statistical Association. Vol. 96, No. 456 (2001).
#'
#' @return A matrix with an attribute df if \code{return.df=TRUE}.
#' 
#' @examples
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#'
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- sequence.lvm(0)
#' set.seed(10)
#' d <- sim(m,n)
#'
#' e.lm <- lm(formula.lvm,data=d)
#' iid2(e.lm)
#' iid2(e.lm,corrected=FALSE)-iid(e.lm)
#' 
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' score(e.lvm)
#' range(iid2(e.lvm, adjust.residuals = FALSE)-iid(e.lvm))
#'
#' #### multiple regression
#' m <- lvm(c(Y1~X1,Y2~X2,Y3~X3))
#' e.lvm <- estimate(m,data=sim(m,1e2))
#' iid2(e.lvm)
#' 
#' @export
`iid2` <-
  function(x, ...) UseMethod("iid2")

## * iid2.lm
#' @rdname iid2
#' @export
iid2.lm <- function(x, data = NULL, adjust.residuals = TRUE, power = 1/2, ...){
    if(!identical(class(x),"lm")){
        wrongClass <- paste(setdiff(class(x),"lm"), collapse = " ")
        stop("iid2 is not available for ",wrongClass," objects \n")
    }
    if(is.null(data)){
        data <- model.matrix(x)
    }else{
        data <- model.matrix(formula(x), data)
    }

### ** get hat values
    XX_m1 <- solve(t(data)%*%data)
    H <- data %*% XX_m1 %*% t(data)    
    if(adjust.residuals){
        epsilon <- residuals(x)/(1 - diag(H))^{power}
    }else{
        epsilon <- residuals(x)
    }

    iid0 <- sweep(data, MARGIN = 1, FUN = "*", STATS = epsilon) %*% XX_m1
    
### ** export
    return(iid0)
}

## * iid2.gls
#' @rdname iid2
#' @export
iid2.gls <- function(x, p = NULL, data = NULL, 
                     adjust.residuals = TRUE, power = 1/2,
                     return.df = TRUE, ...){

### ** normalize argument
    if(is.null(data)){
        data <- getData(x)
    }
    if(is.null(p)){
        p <- coef(x)
    }

### ** compute the score
    e.score <- score2(x, p = p, data = data,
                      adjust.residuals = adjust.residuals, power = power,
                      indiv = TRUE, return.vcov.param = TRUE, ...)
    vcov.param <- attr(e.score,"vcov.param")
    
### ** compute the iid
    iid0 <- e.score %*% vcov.param

### ** compute the degrees of freedom
    if(return.df){
        df.adj <- calcDDF(x, p = p, iid0 = iid0,
                          adjust.residuals = adjust.residuals, power = power, ...)
        attr(iid0,"df") <- df.adj
    }
    
### ** export
    colnames(iid0) <- colnames(e.score)
    
    return(iid0)
}

## * iid2.lme
#' @rdname iid2
#' @export
iid2.lme <- function(x, p = NULL, data = NULL,
                     adjust.residuals = TRUE, power = 1/2,
                     return.df = TRUE, ...){

    if(is.null(data)){
        data <- getData(x)
    }
    if(is.null(p)){
        p <- fixef(x)
    }
    
### ** compute the score
    e.score <- score2(x, p = p, data = data,
                      adjust.residuals = adjust.residuals, power = power,
                      indiv = TRUE, return.vcov.param = TRUE, ...)
    vcov.param <- attr(e.score,"vcov.param")
    
### ** compute the iid
    iid0 <- e.score %*% vcov.param

### ** compute the degrees of freedom
    if(return.df){
        df.adj <- calcDDF(x, p = p, iid0 = iid0,
                          adjust.residuals = adjust.residuals, power = power, ...)
        attr(iid0,"df") <- df.adj
    }
    
### ** export
    colnames(iid0) <- colnames(e.score)
    return(iid0)
}


## * iid2.lvmfit
#' @rdname iid2
#' @export
iid2.lvmfit <- function(x, p = NULL, data = NULL, 
                        adjust.residuals = TRUE, Dmethod = FALSE, power = 1/2, 
                        return.df = TRUE, check.score = TRUE, ...){

    if(is.null(data)){
        data <- model.frame(x)
    }else{
        data <- as.data.frame(data)[,vars(x),drop=FALSE]
    }
    if(is.null(p)){
        p <- coef(x)
    }   
    name.param <- names(p)
    
### ** compute the score
    if(check.score){
        browser()
        S1 <- score2(x, p = p, data = data,
                     adjust.residuals = FALSE, Dmethod = "simple",
                     indiv = TRUE)
        S2 <- score(x, p = p, indiv = TRUE)
        if(max(abs(S1-S2))>1e-10){
            stop("score2 does not match score \n",
                 "report that to the maintainer of the package")
        }
    }
    e.score <- score2(x, p = p, data = data,
                      adjust.residuals = adjust.residuals, power = power,
                      indiv = TRUE, return.vcov.param = TRUE, ...)
    vcov.param <- attr(e.score,"vcov.param")

### ** compute the iiid
    iid0 <- e.score %*% vcov.param

### ** degrees of freedom
    if(return.df){
        browser()
        df.adj <- .calcDF(x = x, p = p, iid0 = iid0)
        attr(iid0,"df") <- df.adj
    }
    
### ** export
    colnames(iid0) <- colnames(e.score)
    return(iid0)
}

##----------------------------------------------------------------------
### iid2.R ends here

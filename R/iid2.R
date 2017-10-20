### iid2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (13:16) 
## Version: 
## last-updated: okt 19 2017 (18:44) 
##           By: Brice Ozenne
##     Update #: 134
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
#' @param ... not used
#' 
#' @references Kauermann G. and Carroll R. J. A note on the efficiency of sandwich covariance matrix estimation. Journal of the American Statistical Association. Vol. 96, No. 456 (2001).
#'
#' @examples
#'
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

## * method iid2.lm
#' @rdname iid2
#' @export
iid2.lm <- function(x, data = NULL, adjust.residuals = TRUE, alpha = 1/2, ...){
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
        epsilon <- residuals(x)/(1 - diag(H)^{alpha})
    }else{
        epsilon <- residuals(x)
    }

    iid0 <- sweep(data, MARGIN = 1, FUN = "*", STATS = epsilon) %*% XX_m1

### ** export
    return(iid0)
}

iid2.lvmfit <- function(x, data = NULL, 
                        adjust.residuals = TRUE, power = 1/2,
                        use.information = FALSE, Dmethod = "Richardson",
                        return.df = TRUE,
                        check.score = TRUE, ...){

    if(!identical(class(x),"lvmfit")){
        wrongClass <- paste(setdiff(class(x),"lvmfit"), collapse = " ")
        stop("iid2 is not available for ",wrongClass," objects \n")
    }

    if(is.null(data)){
        data <- model.frame(x)
    }else{
        stop("not implemented yet!\n")
    }
    
### ** compute the iid
    coef.model <- coef(x)

    if(check.score){
        S1 <- score2(x, param = coef.model, data = data,
                     adjust.residuals = FALSE, indiv = TRUE, return.df = FALSE)
        S2 <- score(x, p = coef.model,indiv = TRUE)
        if(max(abs(S1-S2))>1e-10){
            stop("score2 does not match score \n",
                 "report that to the maintainer of the package")
        }
    }

    e.score <- score2(x, p = coef.model, data = data,
                      adjust.residuals = adjust.residuals, power = power, return.df = return.df)

    if(use.information){
        iI <- vcov(x)
    }else{
        I <- -numDeriv::jacobian(func = function(p){
            score(x, p = p, data = data, indiv = FALSE, ...)        
        },
        x = coef.model,
        method = Dmethod)
        iI <- lava::Inverse(I)
    }

    iid0 <- e.score %*% iI

### ** export
    colnames(iid0) <- colnames(e.score)
    if(return.df){
        attr(iid0,"df") <- attr(e.score,"df")
    }
    return(iid0)
}

##----------------------------------------------------------------------
### iid2.R ends here

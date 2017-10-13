### iid2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (13:16) 
## Version: 
## last-updated: okt 13 2017 (10:26) 
##           By: Brice Ozenne
##     Update #: 76
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
#' 
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' iid2(e.lvm)
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
iid2.lm <- function(x, data = NULL, ...){
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
    epsilon <- residuals(x)/(1 - diag(H))

    iid.corrected <- sweep(data, MARGIN = 1, FUN = "*", STATS = epsilon) %*% XX_m1

### ** export
    return(iid.corrected)
}

iid2.lvmfit <- function(x, data = NULL, ...){

    if(!identical(class(x),"lvmfit")){
        wrongClass <- paste(setdiff(class(x),"lvmfit"), collapse = " ")
        stop("iid2 is not available for ",wrongClass," objects \n")
    }

    if(is.null(data)){
        data <- model.frame(x)
    }else{
        stop("not implemented yet!\n")
    }
    
### ** get empirical moments
    mom.data <- lava:::procdata.lvm(x, data = data)

### ** get theorical moments
    mom.th <- moments(x, p = pars(x), conditional=TRUE, data = data)

### ** compute the score
    e.score <- score2(x, mu = mom.th$xi , Omega =  mom.th$C, data = data)

}

##----------------------------------------------------------------------
### iid2.R ends here

### sCorrect-vcov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 11 2019 (13:55) 
## Version: 
## Last-Updated: jan 27 2020 (13:36) 
##           By: Brice Ozenne
##     Update #: 42
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * vcov2 (documentation)
#' @title  Extract the Variance Covariance Matrix of the Model Parameters After Small Sample Correction
#' @description  Extract the variance covariance matrix of the model parameters from \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} objects.
#' Similar to \code{vcov} but with small sample correction (if any).
#' @name vcov2
#'
#' @param object a \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} object.
#' @param as.lava [logical] Should the order and the name of the coefficient be the same as those obtained using coef with type = -1.
#' Only relevant for \code{lvmfit} objects.
#' 
#' @details If argument \code{p} or \code{data} is not null, then the small sample size correction is recomputed to correct the influence function.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @return A matrix.
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
#' d <- lava::sim(m,n)
#'
#' ## linear model
#' e.lm <- lm(formula.lvm,data=d)
#' vcov.tempo <- vcov2(e.lm)
#' vcov.tempo[rownames(vcov(e.lm)),colnames(vcov(e.lm))]/vcov(e.lm)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' vcov.tempo <- vcov2(e.lvm)
#' vcov.tempo/vcov(e.lvm)
#'
#' @concept small sample inference
#' @export
`vcov2` <-
    function(object, as.lava) UseMethod("vcov2")

## * vcov2.sCorrect
#' @rdname vcov2
#' @export
vcov2.sCorrect <- function(object, as.lava = TRUE){
    out <- .info2vcov(object$sCorrect$information)
    if(as.lava == FALSE){ 
        out <- out[names(object$sCorrect$skeleton$originalLink2param),names(object$sCorrect$skeleton$originalLink2param),drop=FALSE]
        dimames(out) <- list(as.character(object$sCorrect$skeleton$originalLink2param),as.character(object$sCorrect$skeleton$originalLink2param))
    }
    return(out)
}

## * .info2vcov (helper)
#' @title Inverse the Information Matrix
#' @description Compute the inverse of the information matrix.
#' @name vcov2-internal
#'
#' @param attr.info [logical] should the information matrix be returned as an attribute?
#' @param ... arguments passed to .information2
#' 
#' @keywords internal
.info2vcov <- function(information, attr.info = FALSE){
    vcov <- try(chol2inv(chol(information)), silent = TRUE)
    if(inherits(vcov, "try-error")){
        vcov <- try(solve(information), silent = TRUE)
        if(inherits(vcov, "try-error")){ ## try by block
            cat("Singular information matrix: try to inverse it by block \n")
            information.N0 <- abs(information)>1e-10
            remaining.var <- colnames(information)
            vcov <- matrix(0, nrow = NROW(information), ncol = NCOL(information),
                           dimnames = dimnames(information))
            while(length(remaining.var)>0){
                current.set <- remaining.var[1]
                new.set <- unique(unlist(apply(information.N0[current.set,,drop=FALSE],1,function(iRow){list(names(iRow[iRow==1]))})))
                while(length(current.set)<length(new.set)){
                    current.set <- new.set
                    new.set <- unique(unlist(apply(information.N0[current.set,,drop=FALSE],1,function(iRow){list(names(iRow[iRow==1]))})))
                }
                if(length(new.set)>0){
                    iTry <- try(solve(information[current.set,current.set]), silent = TRUE)
                    if(inherits(iTry,"try-error")){
                        vcov[current.set,current.set] <- NA
                    }else{
                        vcov[current.set,current.set] <- iTry
                    }
                }
                remaining.var <- setdiff(remaining.var, current.set)
            }
        }
    }
    if(attr.info){
        attr(vcov,"information") <- information
    }
    if(!inherits(vcov, "try-error")){
        dimnames(vcov) <- dimnames(information)
    }
    return(vcov)
}


######################################################################
### sCorrect-vcov2.R ends here

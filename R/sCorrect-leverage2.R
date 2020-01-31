### leverage2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: feb 19 2018 (17:58) 
## Version: 
## Last-Updated: jan 31 2020 (10:47) 
##           By: Brice Ozenne
##     Update #: 117
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * documentation - leverage2
#' @title Extract Leverage Values
#' @description Extract leverage values from a Gaussian linear model. 
#' @name leverage2
#' 
#' @param object a \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} object.
#' @param format [character] Use \code{"wide"} to return the residuals in the wide format (one row relative to each sample).
#' Otherwise use \code{"long"} to return the residuals in the long format.
#'
#' @details The leverage are defined as the partial derivative of the fitted values with respect to the observations.
#' \deqn{
#' leverage_i = \frac{\partial \hat{Y}_i}{\partial Y_i}
#' }
#' See Wei et al. (1998). \cr \cr
#' 
#' If argument \code{p} or \code{data} is not null, then the small sample size correction is recomputed to correct the residuals.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#' 
#' @return a matrix containing the leverage relative to each sample (in rows)
#' and each endogenous variable (in column).
#'
#' @references Bo-Cheng Wei et al., Generalized Leverage and its applications (1998), Scandinavian Journal of Statistics 25:1:25-37.
#' 
#' @examples
#' ## simulate data
#' set.seed(10)
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#' d <- lava::sim(m,20, latent = FALSE)
#'
#' ## standard linear model
#' e.lm <- lm(Y1~Y2, data = d)
#'
#' sCorrect(e.lm) <- TRUE
#' range(as.double(leverage2(e.lm)) - influence(e.lm)$hat)
#'
#' ## latent variable model
#' e.lvm <- estimate(m, data = d)
#' sCorrect(e.lvm) <- TRUE
#' leverage2(e.lvm)
#' 
#' @concept small sample inference
#' @export
`leverage2` <-
    function(object, format) UseMethod("leverage2")

## * leverage2.sCorrect
#' @rdname leverage2
#' @export
leverage2.sCorrect <- function(object, format = "wide"){

    format <- match.arg(format, choices = c("long","wide"))

    if(format == "wide"){
        return(object$sCorrect$leverage)
    }else if(format == "long"){
        endogenous <- colnames(object$sCorrect$leverage)
        n.endogenous <- length(endogenous)
        
        outW <- data.frame(cluster = 1:NROW(object$sCorrect$leverage), object$sCorrect$leverage)
        outL <- stats::na.omit(stats::reshape(outW,
                                              idvar = "id",
                                              direction = "long",
                                              varying = list(endogenous),
                                              timevar = "endogenous",
                                              v.names = "leverage"))
        rownames(outL) <- NULL
        outL$endogenous <- factor(outL$endogenous, levels = 1:n.endogenous, labels = endogenous)
        reorder <- match(interaction(original.order$XXclusterXX,original.order$XXendogenousXX),
                         interaction(outL$cluster,outL$endogenous))
        return(outL[reorder,])
    }

}

## * .leverage2
.leverage2 <- function(Omega, epsilon, dmu, dOmega, vcov.param,
                       name.pattern, missing.pattern, unique.pattern,
                       endogenous, n.endogenous, param, param.mean, param.hybrid, n.cluster){

    n.pattern <- NROW(unique.pattern)
    n.param <- length(param)
    n.param.mean <- length(param.mean)
    n.param.hybrid <- length(param.hybrid)
    leverage <- matrix(NA, nrow = n.cluster, ncol = n.endogenous,
                       dimnames = list(NULL, endogenous))
    if(is.null(vcov.param)){
        leverage[] <- 0
        return(leverage)
    }

    scoreY <- array(NA, dim = c(n.cluster, n.endogenous, n.param.mean),
                    dimnames = list(NULL, endogenous, param.mean))
    
    for(iP in 1:n.pattern){ ## iP <- 1 
        iIndex <- missing.pattern[[iP]]
        iY <- which(unique.pattern[iP,]==1)
        
        iOmega <- Omega[iY,iY,drop=FALSE]
        iOmegaM1 <- chol2inv(chol(iOmega))
        iOmegaM1.epsilon <- epsilon[iIndex,iY,drop=FALSE] %*% iOmegaM1
            
        ## derivative of the score regarding Y
        for(iParam in param.mean){
            
            if(iParam %in% param.hybrid){
                iOmegaM1.epsilon.dOmega.iOmegaM1 <- iOmegaM1.epsilon %*% dOmega[[iParam]][iY,iY,drop=FALSE] %*% iOmegaM1
            }else{
                iOmegaM1.epsilon.dOmega.iOmegaM1 <- 0
            }

            if(length(iY)>1){
                scoreY[iIndex,iY,iParam] <- t(dmu[iParam,iY,iIndex]) %*% iOmegaM1 + 2 * iOmegaM1.epsilon.dOmega.iOmegaM1
            }else{
                scoreY[iIndex,iY,iParam] <- dmu[iParam,iY,iIndex] * iOmegaM1[1,1] + 2 * iOmegaM1.epsilon.dOmega.iOmegaM1
            }
        }

        ## leverage
        for(iiY in iY){ ## iiY <- iY[2]
            if(n.param.mean==1){                
                leverage[iIndex,iiY] <- dmu[,iiY,iIndex] * vcov.param[param.mean,param.mean] * scoreY[iIndex,iiY,]
            }else{
                leverage[iIndex,iiY] <- rowSums((t(dmu[,iiY,iIndex]) %*% vcov.param[param.mean,param.mean,drop=FALSE]) * scoreY[iIndex,iiY,])
            }
            ## diag( t(dmu[,iiY,iIndex]) %*% vcov.param[param.mean,param.mean,drop=FALSE] %*% t(scoreY[iIndex,iiY,]) )
        }
    }        

    return(leverage)            
}

##----------------------------------------------------------------------
### leverage2.R ends here

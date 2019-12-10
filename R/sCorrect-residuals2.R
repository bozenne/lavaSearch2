### sCorrect-residuals2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (11:17) 
## Version: 
## Last-Updated: dec 10 2019 (16:07) 
##           By: Brice Ozenne
##     Update #: 36
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
#' @title Extract Residuals After Small Sample Correction.
#' @description Extract residuals from \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} objects.
#' Similar to \code{residuals} but with small sample correction.
#' @name residuals2
#' 
#' @param object a \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} object.
#' @param type [character] the type of residual to extract:
#' \code{"response"} for raw residuals,
#' \code{"studentized"} for studentized residuals,
#' \code{"normalized"} for normalized residuals.
#' @param param [named numeric vector] the fitted parameters.
#' @param data [data.frame] the data set.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @details If argument \code{p} or \code{data} is not null, then the small sample size correction is recomputed to correct the residuals. \cr
#'
#' The raw residuals are defined by  observation minus the fitted value:
#' \deqn{
#' \varepsilon = (Y_1 - \mu_1, ..., Y_m - \mu_m)
#' }
#' The studentized residuals divided the raw residuals relative to each endogenous variable by the modeled variance of the endogenous variable.
#' \deqn{
#' \varepsilon_{stud} =(\frac{Y_1 - \mu_1}{\sigma_1}, ..., \frac{Y_m - \mu_m}{\sigma_m})
#' }
#' The normalized residuals multiply the raw residuals by the inverse of the square root of the modeled residual variance covariance matrix.
#' \deqn{
#' \varepsilon_{norm} = \varepsilon \Omega^{-1/2}
#' }
#' @return a matrix containing the residuals relative to each sample (in rows)
#' and each endogenous variable (in column).
#' @concept small sample inference
#' @export
`residuals2` <-
    function(object, param, data, type, format) UseMethod("residuals2")

## * Examples
#' @rdname residuals2
#' @examples
#' set.seed(10)
#' n <- 101
#'
#' #### linear regression ####
#' Y1 <- rnorm(n, mean = 0)
#' Y2 <- rnorm(n, mean = 0.3)
#' Id <- findInterval(runif(n), seq(0.1,1,0.1))
#' data.df <- rbind(data.frame(Y=Y1,G="1",Id = Id),
#'            data.frame(Y=Y2,G="2",Id = Id)
#'            )
#' m.lm <- lm(Y ~ G, data = data.df)
#' residuals2(m.lm)
#'
#' #### mixed models ####
#' library(nlme)
#' m.gls <- gls(Y ~ G, weights = varIdent(form = ~ 1|Id), data = data.df)
#' coef2(m.gls)
#' m.lme <- lme(Y ~ G, random = ~ 1|Id, data = data.df)
#' coef2(m.lme)
#' 
#' #### latent variable models ####
#' library(lava)
#' e.lvm <- estimate(lvm(Y ~ G), data = data.df)
#' coef2(e.lvm)
#'
#' ## simulate data
#' set.seed(10)
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#' d <- lava::sim(m,20, latent = FALSE)
#'
#' ## standard linear model
#' e.lm <- lm(Y1~Y2, data = d)
#' sCorrect(e.lm) <- TRUE
#' 
#' sigma(e.lm)^2
#' mean(residuals(e.lm)^2)
#' mean(residuals2(e.lm)^2)
#' 
#' ## latent variable model
#' e.lvm <- estimate(m, data = d)
#' sCorrect(e.lvm) <- TRUE
#' mean(residuals2(e.lvm)^2)
#'

## * residuals2.lm
#' @rdname residuals2
#' @export
residuals2.lm <- function(object, param = NULL, data = NULL, type = "response", format = "wide"){

    format <- match.arg(format, choices = c("long","wide"))
    
    if(inherits(object,"sCorrect") && is.null(param) && is.null(data)){
        cM <- object$sCorrect
    }else{
        cM <- conditionalMoment(object, param = param, data = data,
                                initialize = !is.null(data) || !inherits(object,"sCorrect"),
                                first.order = FALSE, second.order = FALSE, usefit = TRUE)
    }

    residuals <- .normalizeResiduals(residuals = cM$moment$residuals,
                                     Omega = cM$moment$Omega,
                                     type = type,
                                     missing.pattern = cM$moment$missing.pattern,
                                     Omega.missing.pattern = cM$moment$Omega.missing.pattern)
    if(format == "wide"){
        return(residuals)
    }else if(format == "long"){
        endogenous <- colnames(residuals)
        n.endogenous <- length(endogenous)
        
        residualsW <- data.frame(cluster = 1:NROW(residuals), residuals)
        residualsL <- stats::na.omit(stats::reshape(residualsW,
                                                    idvar = "id",
                                                    direction = "long",
                                                    varying = list(endogenous),
                                                    timevar = "endogenous",
                                                    v.names = "residual"))
        rownames(residualsL) <- NULL
        residualsL$endogenous <- factor(residualsL$endogenous, levels = 1:n.endogenous, labels = endogenous)
        reorder <- match(interaction(cM$original.order$XXclusterXX,cM$original.order$XXendogenousXX),
                         interaction(residualsL$cluster,residualsL$endogenous))
        return(residualsL[reorder,])
    }
}

## * residuals2.gls
#' @rdname residuals2
#' @export
residuals2.gls <- residuals2.lm

## * residuals2.lme
#' @rdname residuals2
#' @export
residuals2.lme <- residuals2.lm

## * residuals2.lvmfit
#' @rdname residuals2
#' @export
residuals2.lvmfit <- residuals2.lm

## * .normalizeResiduals
.normalizeResiduals <- function(residuals, Omega, type,
                                missing.pattern, Omega.missing.pattern){
    type <- match.arg(type, choices = c("response","studentized","normalized"), several.ok = FALSE)

    if(type=="studentized"){
        residuals <- sweep(residuals,
                           STATS = sqrt(diag(Omega)),
                           FUN = "/",
                           MARGIN = 2)
        ## object$sCorrect$residuals/residuals
    }else if(type=="normalized"){
        name.endogenous <- colnames(residuals)
        if(any(is.na(residuals))==FALSE){
            residuals <- residuals %*% matrixPower(Omega, symmetric = TRUE, power = -1/2)
        }else{
            iOmegaHalf.missing.pattern <- lapply(Omega.missing.pattern,matrixPower,symmetric = TRUE, power = -1/2)
            for(iC in 1:NROW(residuals)){ ## iC <- 1
                iY <- which(!is.na(residuals[iC,]))
                residuals[iC,iY] <- residuals[iC,iY] %*% iOmegaHalf.missing.pattern[[missing.pattern[iC]]]
            }
        }
        colnames(residuals) <- name.endogenous
    }

    return(residuals)
}



######################################################################
### sCorrect-residuals2.R ends here

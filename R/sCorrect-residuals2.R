### sCorrect-residuals2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (11:17) 
## Version: 
## Last-Updated: nov 18 2019 (11:22) 
##           By: Brice Ozenne
##     Update #: 5
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
    function(object, param, data, type) UseMethod("residuals2")

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
residuals2.lm <- function(object, param = NULL, data = NULL, type = "response"){
    
    if(is.null(param) && is.null(data)){ ## default value
        Omega <- getVarCov2(object)
        
        residuals <- cbind(stats::residuals(object))
        colnames(residuals) <- colnames(Omega)
    }else if(is.null(param)){ ## new dataset
        response.var <- endogenous(object)
        Omega <- getVarCov2(object)
                                       
        residuals <- cbind(stats::predict(object, newdata = data) - data[[response.var]])
        colnames(residuals) <- colnames(Omega)
    }else{ ## new parameter values (and possibly new dataset)
        stop("Not yet implemented \n")
    }
    
    residuals <- .normalizeResiduals(residuals = residuals,
                                     Omega = Omega,
                                     type = type)
    
    return(residuals)
}

## * residuals2.gls
#' @rdname residuals2
#' @export
residuals2.gls <- function(object, param = NULL, data = NULL, type = "response", ...){
    browser()
    if(is.null(param) && is.null(data)){ ## default value
        Omega <- getVarCov2(object)
        
        residuals <- stats::residuals(object)
        info <- getGroups2(object, ...)

        M.residuals <- matrix(NA, nrow = info$n.cluster, ncol = info$n.endogenous,
                              dimnames = list(NULL,info$name.endogenous))
        M.residuals[info$index.vec2mat] <- residuals
    }else if(is.null(param)){ ## new dataset
        response.var <- endogenous(object)
        Omega <- getVarCov2(object)
                                       
        residuals <- cbind(stats::predict(object, newdata = data) - data[[response.var]])
        colnames(residuals) <- colnames(Omega)
    }else{ ## new parameter values (and possibly new dataset)
        stop("Not yet implemented \n")
    }
    
    residuals <- .normalizeResiduals(residuals = residuals,
                                     Omega = Omega,
                                     type = type)
    
    return(residuals)
}

## * residuals2.lme
#' @rdname residuals2
#' @export
residuals2.lme <- residuals2.gls

## * residuals2.lvmfit
#' @rdname residuals2
#' @export
residuals2.lvmfit <- function(object, param = NULL, data = NULL, type = "response"){

    res <- residuals(object, data = data, p = param)
    return(res)
}

## * residuals2.sCorrect
#' @rdname residuals2
#' @export
residuals2.sCorrect <- function(object, param = NULL, data = NULL, type = "response"){


    if(!is.null(param) || !is.null(data)){
        args <- object$sCorrect$args
        args$df <- FALSE
        object$sCorrect <- do.call(sCorrect,
                                   args = c(list(object, param = param, data = data),
                                            args))
    }
    
    residuals <- .normalizeResiduals(residuals = object$sCorrect$residuals,
                                     Omega = object$sCorrect$Omega,
                                     type = type)
    
    return(residuals)
}

## * .normalizeResiduals
.normalizeResiduals <- function(residuals, Omega, type){
    type <- match.arg(type, choices = c("response","studentized","normalized"), several.ok = FALSE)

    if(type=="studentized"){
        residuals <- sweep(residuals,
                           STATS = sqrt(diag(Omega)),
                           FUN = "/",
                           MARGIN = 2)
        ## object$sCorrect$residuals/residuals
    }else if(type=="normalized"){
        name.endogenous <- colnames(residuals)
        residuals <- residuals %*% matrixPower(Omega, symmetric = TRUE, power = -1/2)
        colnames(residuals) <- name.endogenous
        ## object$sCorrect$residuals/residuals
        ## var(residuals)        
    }

    return(residuals)
}



######################################################################
### sCorrect-residuals2.R ends here

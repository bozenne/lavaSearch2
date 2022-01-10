### Scorrect-confint2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: Jan  4 2022 (10:59) 
## Version: 
## Last-Updated: Jan 10 2022 (15:00) 
##           By: Brice Ozenne
##     Update #: 67
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * Documentation
#' @title Confidence Intervals With Small Sample Correction
#' @description Extract confidence intervals of the coefficients from a latent variable model.
#' Similar to \code{lava::confint} but with small sample correction.
#' @name confint2
#'
#' @param object a \code{lvmfit} or \code{lvmfit2} object (i.e. output of \code{lava::estimate} or \code{lavaSearch2::estimate2}).
#' @param robust [logical] should robust standard errors be used instead of the model based standard errors? Should be \code{TRUE} if argument cluster is not \code{NULL}.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param as.lava [logical] when \code{TRUE} uses the same names as when using \code{stats::coef}.
#' @param transform [function] transformation to be applied.
#' @param conf.level [numeric, 0-1] level of the confidence intervals.
#' @param ssc [character] method used to correct the small sample bias of the variance coefficients: no correction (code{"none"}/\code{FALSE}/\code{NA}),
#' correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
#' Only relevant when using a \code{lvmfit} object. 
#' @param df [character] method used to estimate the degree of freedoms of the Wald statistic: Satterthwaite \code{"satterthwaite"}. 
#' Otherwise (\code{"none"}/code{FALSE}/code{NA}) the degree of freedoms are set to \code{Inf}.
#' Only relevant when using a \code{lvmfit} object. 
#' @param ... additional argument passed to \code{estimate2} when using a \code{lvmfit} object. 
#'
#' @details When argument object is a \code{lvmfit} object, the method first calls \code{estimate2} and then extract the confidence intervals.
#' 
#' @return A data.frame with a row per coefficient.
#'
#' @concept extractor
#' @keywords smallSampleCorrection
#' @export
`confint2` <- function(object, robust, cluster, transform,
                       as.lava, conf.level, ...) UseMethod("confint2")

## * Examples
#' @rdname confint2
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' dW <- sampleRepeated(10, format = "wide")
#' set.seed(10)
#' dL <- sampleRepeated(10, format = "long")
#' dL$time2 <- paste0("visit",dL$time)
#' 
#' #### latent variable models ####
#' e.lvm <- estimate(lvm(c(Y1,Y2,Y3) ~ 1*eta + X1, eta ~ Z1), data = dW)
#' confint(e.lvm)
#' confint2(e.lvm)
#' confint2(e.lvm, as.lava = FALSE)

## * confint2.lvmfit
#' @rdname confint2
confint2.lvmfit <- function(object, robust = FALSE, cluster = NULL,
                            as.lava = TRUE, conf.level = 0.95, transform = NULL,
                            ssc = lava.options()$ssc, df = lava.options()$df, ...){

    return(confint(estimate2(object, ssc = ssc, df = df, dVcov.robust = robust, ...), robust = robust, cluster = cluster, as.lava = as.lava, conf.level = conf.level, transform = transform))

}

## * confint2.lvmfit2
#' @rdname confint2
confint2.lvmfit2 <- function(object, robust = FALSE, cluster = NULL,
                            as.lava = TRUE, conf.level = 0.95, transform = NULL, ...){

    dots <- list(...)
    if(length(dots)>0){
        warning("Argument(s) \'",paste(names(dots),collapse="\' \'"),"\' not used by ",match.call()[1],". \n")
    }

    ## ** new model parameters
    param <- coef(object, as.lava = as.lava)
    name.param <- names(param)
    n.param <- length(name.param)

    ## ** new Wald test
    type <- object$sCorrect$skeleton$type
    type <- type[!is.na(type$originalLink),]
    null <- setNames(rep(0, n.param),name.param)
    if(any(type$detail %in% c("Sigma_var","Psi_var"))){
        param.var <- type[type$detail %in% c("Sigma_var","Psi_var"),"param"]
        null[param.var] <- NA
    }
    table.all <- compare2(object,
                          linfct = name.param,
                          rhs = null,
                          robust = robust,
                          cluster = NULL,
                          F.test = FALSE,
                          as.lava = FALSE)
    tableS.all <- summary(table.all, test = multcomp::adjusted("none"), transform = transform, conf.level = conf.level)$table2
    rownames(tableS.all) <- name.param
    
    if(as.lava){
        out <- tableS.all[names(object$sCorrect$skeleton$originalLink2param),c("lower","upper"),drop=FALSE]
        rownames(out) <- as.character(object$sCorrect$skeleton$originalLink2param)
        colnames(out) <- c(paste0(100*(1-conf.level)/2," %"),paste0(100*(1-(1-conf.level)/2)," %"))
        repturn(out)
    }else{
        return(tableS.all)
    }

}

## *  confint.lvmfit
#' @rdname confint2
#' @export
confint.lvmfit2 <- confint2.lvmfit2

##----------------------------------------------------------------------
### Scorrect-confint2.R ends here

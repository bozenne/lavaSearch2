### sCorrect_init.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  2 2019 (11:52) 
## Version: 
## Last-Updated: dec 12 2019 (14:40) 
##           By: Brice Ozenne
##     Update #: 242
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .convert_sCorrect (generate object)
## ** .convert_sCorrect<-
`.convert_sCorrect<-` <-
    function(object, value) UseMethod("convert_sCorrect<-")

## ** .convert_sCorrect<-.default
`convert_sCorrect<-.default` <- function(object, ..., value){

    ## *** preliminary tests
    if(!inherits(object,"lm") && !inherits(object,"gls") && !inherits(object,"lme") && !inherits(object,"lvmfit")){
        stop("Cannot only convert object of class lm/gls/lme/lvmfit to sCorrect \n")
    }

    ## lvm
    if("multigroupfit" %in% class(object)){
        stop("sCorrect cannot handle multigroup models \n")
    }
    if(inherits(object,"lvmfit") && length(object$model$attributes$ordinal)>0){
        name.t <- names(object$model$attributes$ordinal)
        stop("sCorrect does not handle ordinal variables \n",
             "ordinal variable(s): \"",paste(name.t, collapse = "\" \""),"\"\n")
    }
    
    if(inherits(object,"lvmfit") && length(object$model$attributes$transform)>0){
        name.t <- names(object$model$attributes$transform)
        stop("sCorrect does not handle transformed variables \n",
             "transformed variable(s): \"",paste(name.t, collapse = "\" \""),"\"\n")
    }

    ## *** initialize object
    if(is.null(object$sCorrect) || (value == TRUE)){
        object$sCorrect <- .init_sCorrect(object, ...)
    }else{
        stop("Existing small sample correction \n",
             "Set argument \'value\' to TRUE to overwrite it \n")
    }

    ## *** output
    class(object) <- append("sCorrect",class(object))
    return(object)
}

## ** .init_sCorrect
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' dW <- sampleRepeated(10, format = "wide")
#' set.seed(10)
#' dL <- sampleRepeated(10, format = "long")
#' dL$time2 <- paste0("visit",dL$time)
#' 
#' #### linear model ####
#' e.lm <- lm(Y1~X1, data = dW)
#' .init_sCorrect(e.lm)
#'
#' #### gls model ####
#' e.gls1 <- gls(Y1~X1, data = dW)
#' .getGroups2(e.gls1, data = dW)
#' 
#' e.gls2 <- gls(Y~X1, correlation = corCompSymm(form=~1|id), data = dL)
#' .getGroups2(e.gls2, data = dL)
#'
#' e.gls3 <- gls(Y~X1, correlation = corCompSymm(form=~time|id), data = dL)
#' .getGroups2(e.gls3, data = dL)
#'
#' e.gls4 <- gls(Y~X1, weight = varIdent(form=~1|time2), data = dL)
#' .getGroups2(e.gls4, data = dL)
#'
#' e.gls5 <- gls(Y~X1, weight = varIdent(form=~1|time2),
#'               correlation = corSymm(form=~time|id), data = dL)
#' .getGroups2(e.gls5, data = dL)
#'
#' #### lme model ####
#' e.lme <- lme(Y~X1, random=~1|id, data = dL)
#' .getGroups2(e.lme, data = dL)
#' 

.init_sCorrect <- function(object, param = NULL, data = NULL,
                           hessian = TRUE, dInformation = TRUE, ...){
    if(TRUE){cat(".init_sCorrect \n")}

    object$sCorrect <- list(ssc = FALSE, df = FALSE)

    ## ** type of variables
    object$sCorrect$endogenous <- endogenous(object)
    object$sCorrect$exogenous <- exogenous(object)
    object$sCorrect$latent <- latent(object)

    ## ** model coefficients
    if(is.null(param)){
        object$sCorrect$coef <- coef2(object, ssc = FALSE, labels = 1)
    }else{
        object$sCorrect$coef <- param
    }
    
    ## ** data
    if(is.null(data)){
        object$sCorrect$data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                                            envir = parent.env(environment()), rm.na = TRUE)
    }else{
        object$sCorrect$data <- data
    }

    ## ** design matrix
    object$sCorrect$X <- .getDesign(object, data = object$sCorrect$data, add.Y = TRUE)

    ## ** clusters
    object$sCorrect$cluster <- .getGroups2(object, data = object$sCorrect$data, endogenous = object$sCorrect$endogenous)

    ## ** sample size
    object$sCorrect$n.corrected <- object$sCorrect$cluster$n.cluster

    ## ** moments
    object$sCorrect <- c(object$sCorrect,
                         conditionalMoment(object,
                                           initialize = TRUE, first.order = TRUE, second.order = (hessian || dInformation), usefit = TRUE,
                                           param = object$sCorrect$coef,
                                           X = object$sCorrect$X,
                                           cluster = object$sCorrect$cluster,
                                           endogenous = object$sCorrect$endogenous,
                                           latent = object$sCorrect$latent))

    ## ** leverage
    object$sCorrect$leverage <- matrix(0, nrow = object$sCorrect$cluster$n.cluster, ncol = length(object$sCorrect$endogenous),
                                       dimnames = list(NULL, object$sCorrect$endogenous))
    object$sCorrect$leverage[is.na(object$sCorrect$moment$residuals)] <- NA

    ## ** score
    object$sCorrect$score <- .score2(dmu = object$sCorrect$dmoment$dmu,
                                     dOmega = object$sCorrect$dmoment$dOmega,                    
                                     epsilon = object$sCorrect$moment$residuals,
                                     OmegaM1 = object$sCorrect$moment$iOmega.missing.pattern,
                                     missing.pattern = object$sCorrect$missing$pattern,
                                     unique.pattern = object$sCorrect$missing$unique.pattern,
                                     name.pattern = object$sCorrect$missing$name.pattern,
                                     name.param = object$sCorrect$skeleton$Uparam,
                                     name.meanparam = object$sCorrect$skeleton$Uparam.mean,
                                     name.varparam = object$sCorrect$skeleton$Uparam.var,
                                     n.cluster = object$sCorrect$cluster$n.cluster)
    
    ## ** information matrix
    object$sCorrect$information <- .information2(dmu = object$sCorrect$dmoment$dmu,
                                                 dOmega = object$sCorrect$dmoment$dOmega,
                                                 OmegaM1 = object$sCorrect$moment$iOmega.missing.pattern,
                                                 missing.pattern = object$sCorrect$missing$pattern,
                                                 unique.pattern = object$sCorrect$missing$unique.pattern,
                                                 name.pattern = object$sCorrect$missing$name.pattern,
                                                 grid.mean = object$sCorrect$skeleton$grid.dmoment$mean, 
                                                 grid.var = object$sCorrect$skeleton$grid.dmoment$var, 
                                                 name.param = object$sCorrect$skeleton$Uparam,
                                                 leverage = object$sCorrect$leverage,
                                                 n.cluster = object$sCorrect$cluster$n.cluster)
    object$sCorrect$vcov.param <- .info2vcov(object$sCorrect$information, attr.info = FALSE)
    
    ## ** hessian
    if(hessian){
        object$sCorrect$hessian <- .hessian2(dmu = object$sCorrect$dmoment$dmu,
                                             dOmega = object$sCorrect$dmoment$dOmega,
                                             d2mu = object$sCorrect$dmoment$d2mu,
                                             d2Omega = object$sCorrect$dmoment$d2Omega,
                                             epsilon = object$sCorrect$moment$residuals,                                     
                                             OmegaM1 = object$sCorrect$moment$iOmega.missing.pattern,
                                             missing.pattern = object$sCorrect$missing$pattern,
                                             unique.pattern = object$sCorrect$missing$unique.pattern,
                                             name.pattern = object$sCorrect$missing$name.pattern,
                                             grid.mean = object$sCorrect$skeleton$grid.dmoment$mean, 
                                             grid.var = object$sCorrect$skeleton$grid.dmoment$var, 
                                             grid.hybrid = object$sCorrect$skeleton$grid.dmoment$hybrid, 
                                             name.param = object$sCorrect$skeleton$Uparam,
                                             leverage = object$sCorrect$leverage,
                                             n.cluster = object$sCorrect$cluster$n.cluster)
    }

    ## ** dInformation
    if(dInformation){
        browser()
        object$sCorrect$dInformation <- .dInformation2(dmu = object$sCorrect$dmoment$dmu,
                                                       dOmega = object$sCorrect$dmoment$dOmega,
                                                       epsilon = object$sCorrect$moment$residuals,                                     
                                                       OmegaM1 = object$sCorrect$moment$iOmega.missing.pattern,
                                                       missing.pattern = object$sCorrect$missing$pattern,
                                                       unique.pattern = object$sCorrect$missing$unique.pattern,
                                                       name.pattern = object$sCorrect$missing$name.pattern,
                                                       grid.meanparam = object$sCorrect$skeleton$grid$mean, 
                                                       grid.varparam = object$sCorrect$skeleton$grid$var, 
                                                       grid.hybrid = object$sCorrect$skeleton$grid$hybrid, 
                                                       name.param = object$sCorrect$skeleton$Uparam,
                                                       leverage = object$sCorrect$leverage,
                                                       n.cluster = object$sCorrect$cluster$n.cluster)
    }
    
    ## ** output
    return(object)
}


######################################################################
### object-sCorrect.R ends here

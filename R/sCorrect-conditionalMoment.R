### conditionalMoment.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: dec 10 2019 (16:37) 
##           By: Brice Ozenne
##     Update #: 1334
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * conditionalMoment - documentation
#' @title Prepare the Computation of score2
#' @description Compute the conditional mean and variance,
#' and their first and second derivative regarding the model parameters.
#' @name conditionalMoment
#' 
#' @param object,x a latent variable model.
#' @param data [data.frame] data set.
#' @param param [numeric vector] the fitted coefficients.
#' @param cluster [list] information about the cluster and endogenous variables (output of getGroups2).
#' @param first.order [logical] should the terms relative to the first derivative of the conditional moments be considered?
#' @param second.order [logical] should the terms relative to the second derivative of the conditional moments be considered?
#' @param usefit [logical] If TRUE the coefficients estimated by the model are used to pre-compute quantities. 
#' @param ... [internal] only used by the generic method or by the <- methods.
#' 
#' @details For lvmfit objects, there are two levels of pre-computation:
#' \itemize{
#' \item a basic one that do no involve the model coefficient (\code{conditionalMoment.lvm}).
#' \item an advanced one that require the model coefficients (\code{conditionalMoment.lvmfit}). 
#' }
#' 
#' @examples
#' m <- lvm(Y1~eta,Y2~eta,Y3~eta)
#' latent(m) <- ~eta
#'
#' d <- lava::sim(m,1e2)
#' e <- estimate(m, d)
#'
#' ## basic pre-computation
#' res1 <- conditionalMoment(e, data = d,
#'                          first.order = FALSE, second.order = FALSE,
#'                          name.endogenous = endogenous(e),
#'                          name.latent = latent(e), usefit = FALSE)
#' res1$skeleton$Sigma
#' 
#' ## full pre-computation
#' res2 <- conditionalMoment(e, param = coef(e), data = d,
#'                          first.order = FALSE, second.order = FALSE,
#'                          name.endogenous = endogenous(e),
#'                          name.latent = latent(e), usefit = TRUE
#' )
#' res2$value$Sigma
#'
#' @concept small sample inference
#' @concept derivative of the score equation
#' 
#' @keywords internal
#' @export
`conditionalMoment` <-
  function(object, ...) UseMethod("conditionalMoment")

## * conditionalMoment.lm
#' @rdname conditionalMoment
#' @export
conditionalMoment.lm <- function(object,
                                 initialize, first.order, second.order, usefit,
                                 param = NULL,
                                 data = NULL,
                                 X = NULL,
                                 cluster = NULL,
                                 endogenous = NULL,
                                 latent = NULL,
                                 ...){
    if(TRUE){cat("conditionalMoment \n")}

    out <- list()
    
    ## ** get information from object
    if(is.null(endogenous)){
        endogenous <- lava::endogenous(object, format = "wide")
    }
    n.endogenous <- length(endogenous)
    if(is.null(latent)){
        latent <- lava::latent(object)
    }
    if(is.null(X) && !is.null(cluster)){
        stop("Arguments \'X\' and \'cluster\' must either be both specified or both set to NULL \n")
    }
    if(is.null(data) && (is.null(X) || is.null(cluster))){
        data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                            envir = parent.env(environment()), rm.na = TRUE)
    }
    if(is.null(X)){
        X <- .getDesign(object, data = data, add.Y = TRUE)
    }
    if(is.null(cluster)){
        cluster <- .getGroups2(object, data = data, endogenous = endogenous)
    }

    if(any(colnames(X) %in% c("XXvalueXX","XXendogenousXX","XXclusterXX"))){
        stop("\"XXvalueXX\", \"XXendogenousXX\", and \"XXclusterXX\" should not correspond to variable names \n",
             "It is used internally for data manipulation \n")
    }
    if(inherits(object,"lvmfit") || inherits(object,"lvm")){## convert to long format
        X.latent <- matrix(NA, nrow = NROW(X), ncol = length(latent),
                           dimnames = list(NULL, latent))

        X.long <- melt(data.frame(X,X.latent, XXclusterXX = unique(cluster$index.cluster)),
                       id.vars = c("XXclusterXX",setdiff(colnames(X),c(endogenous,latent))),
                       measure.vars = c(endogenous,latent),
                       variable.name = "XXendogenousXX",
                       value.name = "XXvalueXX")

    }else{
        X.long <- X
        names(X.long)[names(X.long) == lava::endogenous(object, format = "long")] <- "XXvalueXX"
        X.long$XXendogenousXX <- NA
        for(iC in 1:cluster$n.cluster){ ## iC <- 10
            X.long$XXendogenousXX[cluster$index.cluster == iC] <- endogenous[na.omit(cluster$index.Omega[[iC]])]
        }
        X.long <- cbind(X.long, XXclusterXX = cluster$index.cluster)
    }
    out$original.order <- X.long[,c("XXclusterXX","XXendogenousXX")]
    X.long <- X.long[order(X.long$XXclusterXX,X.long$XXendogenousXX),]

    ## ** sanity checks
    if(first.order == FALSE && second.order == TRUE){
        stop("Cannot pre-compute second order derivatives without the first order derivatives. \n")
    }
    if(initialize == FALSE && is.null(object$sCorrect$conditionalMoment)){
        stop("Initialization of the conditional moments missing. \n",
             "Consider setting the argument \'initialize\' to TRUE \n")
    }
    
    ## ** initialize
    if(initialize){
        ## *** initialize conditional moments
        out$skeleton <- skeleton(object, X = X.long,
                                 endogenous = endogenous, latent = latent,
                                 n.cluster = cluster$n.cluster,
                                 index.Omega = cluster$index.Omega)

        ## *** initialize partial derivatives of the conditional moments
        if(first.order){
            out$skeleton <- skeletonDtheta(out$skeleton,
                                           X = X.long,
                                           endogenous = endogenous, latent = latent,
                                           n.cluster = cluster$n.cluster,
                                           index.Omega = cluster$index.Omega)
        }
    
        ## *** initialize second order partial derivatives of the conditional moments
        if(second.order){
            out$skeleton <- skeletonDtheta2(moment = out$skeleton,
                                            X = X.long,
                                            endogenous = endogenous, latent = latent,
                                            n.cluster = cluster$n.cluster,
                                            index.Omega = cluster$index.Omega)
        }
    }

    ## ** update according to the value of the model coefficients
    if(usefit){
        if(is.null(param)){
            param <- coef2(object)
        }
        
        ## *** conditional moments
        out$moment <- updateMoment(skeleton = out$skeleton$param,
                                   value = out$skeleton$value,
                                   toUpdate = out$skeleton$toUpdate.moment,
                                   param = param,
                                   endogenous = endogenous,
                                   latent = latent,
                                   n.cluster = cluster$n.cluster)

        ## *** compute residuals
        out$moment$residuals <- matrix(NA, nrow = cluster$n.cluster, ncol = n.endogenous,
                                       dimnames = list(NULL, endogenous))

        for(iEndo in 1:n.endogenous){ ## iEndo <- 1
            iIndexLong <- out$skeleton$obsByEndoInX[[endogenous[iEndo]]]
            iIndexWide <- X.long[iIndexLong,"XXclusterXX"]
            out$moment$residuals[iIndexWide,endogenous[iEndo]] <- X.long[iIndexLong,"XXvalueXX"] - out$moment$mu[iIndexWide,endogenous[iEndo]]
        }

        ## *** identify missing patterns
        missing.pattern <- out$moment$residuals
        missing.pattern[!is.na(missing.pattern)] <- 1
        missing.pattern[is.na(missing.pattern)] <- 0
        Umissing.pattern <- unique(missing.pattern)

        out$moment$missing.pattern <- apply(missing.pattern, MARGIN = 1, FUN = paste0, collapse="")
        out$moment$Umissing.pattern <- Umissing.pattern
        out$moment$Omega.missing.pattern <- lapply(1:NROW(Umissing.pattern), function(iM){ ## iM <- 1
            iIndex <- which(missing.pattern[iM,]==1)
            return(out$moment$Omega[iIndex,iIndex,drop=FALSE])
        })
        names(out$moment$Omega.missing.pattern) <- apply(Umissing.pattern, MARGIN = 1, FUN = paste0, collapse="")
        out$moment$iOmega.missing.pattern <- lapply(out$moment$Omega.missing.pattern, solve)
        
        ## *** first order derivatives
        if(first.order){            
            out$dmoment <- updateDMoment(moment = out$moment,
                                         skeleton = out$skeleton,
                                         toUpdate = out$skeleton$toUpdate.dmoment,
                                         endogenous = endogenous,
                                         latent = latent,
                                         n.cluster = cluster$n.cluster)
        }

        
        ## *** second order derivatives
        if(second.order){
            browser()
            out$d2moment <- skeletonDtheta2(object)
            object$conditionalMoment$d2mu <- out2$d2mu
            object$conditionalMoment$d2Omega <- out2$d2Omega
        }
       
    }

    
    ## ** export
    return(out)
}

## * conditionalMoment.gls
#' @rdname conditionalMoment
#' @export
conditionalMoment.gls <- conditionalMoment.lm

## * conditionalMoment.lme
#' @rdname conditionalMoment
#' @export
conditionalMoment.lme <- conditionalMoment.lm

## * conditionalMoment.lvmfit
#' @rdname conditionalMoment
#' @export
conditionalMoment.lvmfit <- conditionalMoment.lm


### conditionalMoment.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: jan  7 2020 (14:02) 
##           By: Brice Ozenne
##     Update #: 1439
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
    function(object,
             initialize, first.order, second.order, usefit,
             param = NULL,
             data = NULL) UseMethod("conditionalMoment")

## * conditionalMoment.lm
#' @rdname conditionalMoment
#' @export
conditionalMoment.lm <- function(object,
                                 initialize, first.order, second.order, usefit,
                                 param = NULL,
                                 data = NULL){
    if(lava.options()$debug){cat("conditionalMoment \n")}

    ## ** sanity checks
    if(first.order == FALSE && second.order == TRUE){
        stop("Cannot pre-compute second order derivatives without the first order derivatives. \n")
    }
    if(initialize == FALSE && is.null(object$sCorrect)){
        stop("Initialization of the conditional moments missing. \n",
             "Consider setting the argument \'initialize\' to TRUE \n")
    }

    ## ** initialize
    if(initialize){
        out <- list()

        ## *** get information from object
        out$endogenous <- lava::endogenous(object, format = "wide")
        out$latent <- lava::latent(object)

        if(is.null(data)){
            out$data <- extractData(object, design.matrix = FALSE, as.data.frame = TRUE,
                                    envir = parent.env(environment()), rm.na = TRUE)
        }else{
            out$data <- data
        }
        
        out$X <- .getDesign(object, data = out$data, add.Y = TRUE)
        out$cluster <- .getGroups2(object, data = out$data, endogenous = out$endogenous)

        if(any(colnames(out$X) %in% c("XXvalueXX","XXendogenousXX","XXclusterXX"))){
            stop("\"XXvalueXX\", \"XXendogenousXX\", and \"XXclusterXX\" should not correspond to variable names \n",
                 "It is used internally for data manipulation \n")
        }

        ## *** reshape dataset
        if(inherits(object,"lvmfit") || inherits(object,"lvm")){## convert to long format
            X.latent <- matrix(NA, nrow = NROW(out$X), ncol = length(out$latent),
                               dimnames = list(NULL, out$latent))

            X.long <- melt(data.frame(out$X,X.latent, XXclusterXX = unique(out$cluster$index.cluster)),
                           id.vars = c("XXclusterXX",setdiff(colnames(out$X),c(out$endogenous,out$latent))),
                           measure.vars = c(out$endogenous,out$latent),
                           variable.name = "XXendogenousXX",
                           value.name = "XXvalueXX")

            X.wide <- out$X[,out$endogenous,drop=FALSE]
        }else{
            X.long <- out$X
            names(X.long)[names(X.long) == lava::endogenous(object, format = "long")] <- "XXvalueXX"
            X.long$XXendogenousXX <- NA
            for(iC in 1:out$cluster$n.cluster){ ## iC <- 10
                X.long$XXendogenousXX[out$cluster$index.cluster == iC] <- out$endogenous[na.omit(out$cluster$index.Omega[[iC]])]
            }
            X.long <- cbind(X.long, XXclusterXX = out$cluster$index.cluster)

            X.wide <- dcast(X.long[X.long$XXendogenousXX %in% out$endogenous,c("XXclusterXX","XXendogenousXX","XXvalueXX")],
                            value.var = "XXvalueXX",
                            formula = XXclusterXX ~ XXendogenousXX)
            X.wide <- X.wide[order(X.wide$XXclusterXX),out$endogenous,drop=FALSE]
        }
        out$old2new.order <- X.long[,c("XXclusterXX","XXendogenousXX","XXvalueXX")]
        names(out$old2new.order) <- c("XXclusterXX.old","XXendogenousXX.old","XXvalueXX.old")
        X.long <- X.long[order(X.long$XXclusterXX,X.long$XXendogenousXX),]
        out$old2new.order$XXclusterXX.new <- out$old2new.order$XXclusterXX
        out$old2new.order$XXendogenousXX.new <- out$old2new.order$XXendogenousXX
        out$old2new.order$XXvalueXX.new <- out$old2new.order$XXvalueXX

        ## *** identify missing pattern
        pattern <- X.wide
        pattern[!is.na(pattern)] <- 1
        pattern[is.na(pattern)] <- 0
        unique.pattern <- unique(pattern)
        name.pattern <- apply(unique.pattern, MARGIN = 1, FUN = paste0, collapse = "")
        rownames(unique.pattern) <- name.pattern

        out$missing$pattern <- tapply(1:out$cluster$n.cluster,apply(pattern, MARGIN = 1, FUN = paste0, collapse=""),list)
        out$missing$unique.pattern <- unique.pattern
        out$missing$name.pattern <- name.pattern

        ## *** initialize conditional moments
        out$skeleton <- skeleton(object, X = X.long,
                                 endogenous = out$endogenous, latent = out$latent,
                                 n.cluster = out$cluster$n.cluster,
                                 index.Omega = out$cluster$index.Omega)

        ## *** initialize partial derivatives of the conditional moments
        if(first.order){
            out$skeleton <- skeletonDtheta(out$skeleton,
                                           X = X.long,
                                           endogenous = out$endogenous, latent = out$latent,
                                           missing.pattern = out$missing$pattern,
                                           unique.pattern = out$missing$unique.pattern,
                                           name.pattern = out$missing$name.pattern,
                                           n.cluster = out$cluster$n.cluster,
                                           index.Omega = out$cluster$index.Omega)
        }
    
        ## *** initialize second order partial derivatives of the conditional moments
        if(second.order){
            out$skeleton <- skeletonDtheta2(out$skeleton)
        }
        
    }else{
        out <- object$sCorrect[c("endogenous", "latent", "data", "X", "cluster", "old2new.order", "missing", "skeleton")]
    }

    ## ** update according to the value of the model coefficients
    if(usefit){
        if(is.null(param)){
            out$param <- .coef2(object, labels = 1)[out$skeleton$Uparam]
        }else{
            if(any(names(param) %in% out$skeleton$Uparam == FALSE)){
                stop("Incorrect name for the model parameters \n",
                     "Consider using the argument \'label\' set to 1 when calling stats::coef on the LVM \n")
            }
            out$param <- param[out$skeleton$Uparam]
        }

        ## *** conditional moments
        out$moment <- updateMoment(skeleton = out$skeleton$param,
                                   value = out$skeleton$value,
                                   toUpdate = out$skeleton$toUpdate.moment,
                                   param = out$param,
                                   name.pattern = out$missing$name.pattern,
                                   unique.pattern = out$missing$unique.pattern,
                                   endogenous = out$endogenous,
                                   latent = out$latent,
                                   n.cluster = out$cluster$n.cluster)
        
        ## *** first order derivatives
        if(first.order){            
            out$dmoment <- updateDMoment(moment = out$moment,
                                         skeleton = out$skeleton,
                                         param = out$param)
        }

        
        ## *** second order derivatives
        if(second.order){
            out$d2moment <- updateD2Moment(moment = out$moment,
                                           skeleton = out$skeleton,
                                           param = out$param)
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


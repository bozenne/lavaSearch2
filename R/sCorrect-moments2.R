### conditionalMoment.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 27 2017 (16:59) 
## Version: 
## last-updated: Jan  4 2022 (17:22) 
##           By: Brice Ozenne
##     Update #: 1707
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * moments2 - documentation
#' @title Compute Key Quantities of a Latent Variable Model
#' @description Compute conditional mean, conditional variance, their first and second derivative regarding model parameters, as well as various derivatives of the log-likelihood.
#' @name moments2
#' 
#' @param object a latent variable model.
#' @param data [data.frame] dataset if different from the one used to fit the model.
#' @param param [numeric vector] value of the model parameters if different from the estimated ones.
#' @param initialize [logical] Pre-compute quantities dependent on the data but not on the parameters values.
#' @param usefit [logical] Compute key quantities based on the parameter values.
#' @param score [logical] should the score be output?
#' @param information [logical] should the expected information be output?
#' @param hessian [logical] should the hessian be output?
#' @param vcov [logical] should the variance-covariance matrix based on the expected information be output?
#' @param dVcov [logical] should the derivative of the variance-covariance matrix be output?
#' @param dVcov.robust [logical]  should the derivative of the robust variance-covariance matrix be output?
#' @param Psi [matrix]  Average first order bias in the residual variance. Only necessary for computing adjusted residuals.
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
#' res1 <- conditionalMoment(e, data = d, initialize = TRUE,
#'                          first.order = FALSE, second.order = FALSE, usefit = FALSE)
#' res1$skeleton$param$Sigma
#' 
#' ## full pre-computation
#' res2 <- conditionalMoment(e, param = coef(e), data = d, initialize = TRUE,
#'                          first.order = FALSE, second.order = FALSE, usefit = TRUE)
#' res2$value$Sigma
#'
#' @concept small sample inference
#' @concept derivative of the score equation
#' 
#' @keywords internal
#' @export
`moments2` <-
    function(object, data = NULL, param = NULL, Psi = NULL,
             initialize, usefit,
             score, information, hessian, vcov, dVcov, dVcov.robust, residuals, leverage, derivative) UseMethod("moments2")

## * moments2.lvm
#' @rdname moments
#' @export
moments2.lvm <- function(object, param = NULL, data = NULL, Psi = NULL,
                         initialize, usefit,
                         score, information, hessian, vcov, dVcov, dVcov.robust, residuals, leverage,
                         derivative){
    if(lava.options()$debug){cat("moments2 \n")}
    
    ## ** sanity checks
    if(initialize == FALSE && is.null(object$sCorrect)){
        stop("Initialization of the moments missing. \n",
             "Consider setting the argument \'initialize\' to TRUE \n")
    }
    derivative <- match.arg(derivative, choices = c("analytic","numeric"))

    ## ** initialize
    if(is.null(hessian)){
        hessian <- (derivative == "analytic") && dVcov.robust
    }
    if(score || information || hessian || vcov || dVcov || dVcov.robust){
        first.order <- TRUE
    }else{
        first.order <- FALSE
    }
    if(hessian || dVcov || dVcov.robust){
        second.order <- TRUE
    }else{
        second.order <- FALSE
    }

    ## NOTE: the estimation of the leverage depends on the information matrix and vice versa.
    ##       moments2 is a single run function called by estimate2 which take care of iterating until reaching a stable point
    previous.vcov.param <- object$sCorrect$vcov.param

    ## ** skeleton
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
        out$X <- as.matrix(out$data[,lava::manifest(object),drop=FALSE])
        out$cluster <- .getGroups2(object, data = out$data, endogenous = out$endogenous)

        if(any(colnames(out$X) %in% c("XXvalueXX","XXendogenousXX","XXendogenousXX.index","XXclusterXX"))){
            stop("\"XXvalueXX\", \"XXendogenousXX\", \"XXendogenousXX.index\", and \"XXclusterXX\" should not correspond to variable names \n",
                 "It is used internally for data manipulation \n")
        }

        ## *** reshape dataset (convert to long format)
        X.latent <- matrix(NA, nrow = NROW(out$X), ncol = length(out$latent),
                           dimnames = list(NULL, out$latent))

        X.long <- melt(data.frame(out$X,X.latent, XXclusterXX = unique(out$cluster$index.cluster)),
                       id.vars = c("XXclusterXX",setdiff(colnames(out$X),c(out$endogenous,out$latent))),
                       measure.vars = c(out$endogenous,out$latent),
                       variable.name = "XXendogenousXX",
                       value.name = "XXvalueXX")

        X.wide <- out$X[,out$endogenous,drop=FALSE]

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

        out$missing$pattern <- tapply(1:out$cluster$n.cluster,apply(pattern, MARGIN = 1, FUN = paste0, collapse=""),list)[name.pattern]
        out$missing$unique.pattern <- unique.pattern
        out$missing$name.pattern <- name.pattern

        ## *** initialize conditional moments
        out$skeleton <- skeleton(object, X = X.long,
                                 endogenous = out$endogenous, latent = out$latent,
                                 n.cluster = out$cluster$n.cluster,
                                 index.Omega = out$cluster$index.Omega)

        ## *** initialize partial derivatives of the conditional moments
        out$skeleton <- skeletonDtheta(out$skeleton,
                                       X = X.long,
                                       endogenous = out$endogenous, latent = out$latent,
                                       missing.pattern = out$missing$pattern,
                                       unique.pattern = out$missing$unique.pattern,
                                       name.pattern = out$missing$name.pattern,
                                       n.cluster = out$cluster$n.cluster,
                                       index.Omega = out$cluster$index.Omega)

        ## *** initialize second order partial derivatives of the conditional moments
        ## GS <- skeletonDtheta2(out$skeleton)
        out$skeleton <- skeletonDtheta2(out$skeleton)

        ## *** weights
        if(!is.null(object$weights)){
            out$weights <- object$weights[,1]
        }

    }else{
        ## subset to remove existing results
        rm.name <- c("moment","dmoment","d2moment","score","vcov.param","information","hessian","dInformation","dVcov.param","dRvcov.param","leverage","residuals")
        out <- object$sCorrect[setdiff(names(object$sCorrect),rm.name)]
    }

    ## ** update according to the value of the model coefficients
    if(usefit){
        if(is.null(param)){
            param.tempo <- stats::coef(object, type = 2, labels = 1)
            out$param <- setNames(param.tempo[,"Estimate"],rownames(param.tempo))[out$skeleton$Uparam]
            out$name.param <- setNames(out$skeleton$type[!is.na(out$skeleton$type$originalLink),"param"],
                                       out$skeleton$type[!is.na(out$skeleton$type$originalLink),"originalLink"])
            ## out$name.param <- out$name.param[names(stats::coef(object))]
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

        ## *** update residuals
        if(residuals){
            out$residuals <- .adjustResiduals(epsilon = out$skeleton$param$endogenous - out$moment$mu,
                                              Omega = out$moment$Omega, Psi = Psi, ## Note: if Psi is null returns epsilon i.e. no adjustment
                                              name.pattern = out$missing$name.pattern, missing.pattern = out$missing$pattern, unique.pattern = out$missing$unique.pattern,
                                              endogenous = out$endogenous, n.cluster = out$cluster$n.cluster)
        }
        ## mean(out$residuals^2)
        ## out$moment$Omega
        
        ## *** score
        if(score){
            out$score <- .score2(dmu = out$dmoment$dmu,
                                 dOmega = out$dmoment$dOmega,                    
                                 epsilon = out$residuals,
                                 OmegaM1 = out$moment$OmegaM1.missing.pattern,
                                 missing.pattern = out$missing$pattern,
                                 unique.pattern = out$missing$unique.pattern,
                                 name.pattern = out$missing$name.pattern,
                                 name.param = out$skeleton$Uparam,
                                 name.meanparam = out$skeleton$Uparam.mean,
                                 name.varparam = out$skeleton$Uparam.var,
                                 n.cluster = out$cluster$n.cluster,
                                 weights = out$weights)
        }
    
       
        ## *** leverage
        if(leverage){
            out$leverage <- .leverage2(Omega = out$moment$Omega, 
                                       epsilon = out$residuals,
                                       dmu = aperm(abind::abind(out$dmoment$dmu, along = 3), perm = c(3,2,1)),
                                       dOmega = out$dmoment$dOmega,
                                       vcov.param = previous.vcov.param,
                                       name.pattern = out$missing$name.pattern,
                                       missing.pattern = out$missing$pattern,
                                       unique.pattern = out$missing$unique.pattern,
                                       endogenous = out$endogenous,
                                       n.endogenous = length(out$endogenous),
                                       param = out$skeleton$Uparam,
                                       param.mean = out$skeleton$Uparam.mean,
                                       param.hybrid = out$skeleton$Uparam.hybrid,
                                       n.cluster = out$cluster$n.cluster)
        }

        ## *** information matrix
        if(information){
            out$information <- .information2(dmu = out$dmoment$dmu,
                                             dOmega = out$dmoment$dOmega,
                                             OmegaM1 = out$moment$OmegaM1.missing.pattern,
                                             missing.pattern = out$missing$pattern,
                                             unique.pattern = out$missing$unique.pattern,
                                             name.pattern = out$missing$name.pattern,
                                             grid.mean = out$skeleton$grid.dmoment$mean, 
                                             grid.var = out$skeleton$grid.dmoment$var, 
                                             name.param = out$skeleton$Uparam,
                                             leverage = out$leverage,
                                             n.cluster = out$cluster$n.cluster,
                                             weights = out$weights)
        }
        if(vcov){
            out$vcov.param  <- .info2vcov(out$information, attr.info = FALSE)
        }
    
        ## *** hessian
        if(hessian){
            out$hessian <- .hessian2(dmu = out$dmoment$dmu,
                                     dOmega = out$dmoment$dOmega,
                                     d2mu = out$d2moment$d2mu,
                                     d2Omega = out$d2moment$d2Omega,
                                     epsilon = out$residuals,                                     
                                     OmegaM1 = out$moment$OmegaM1.missing.pattern,
                                     missing.pattern = out$missing$pattern,
                                     unique.pattern = out$missing$unique.pattern,
                                     name.pattern = out$missing$name.pattern,
                                     grid.mean = out$skeleton$grid.dmoment$mean, 
                                     grid.var = out$skeleton$grid.dmoment$var, 
                                     grid.hybrid = out$skeleton$grid.dmoment$hybrid, 
                                     name.param = out$skeleton$Uparam,
                                     leverage = out$leverage,
                                     n.cluster = out$cluster$n.cluster,
                                     weights = out$weights)
        }

        ## *** dVcov.param (model based variance, analytic)
        if((derivative == "analytic") && (dVcov || dVcov.robust)){
            out$dInformation <- .dInformation2(dmu = out$dmoment$dmu,
                                               dOmega = out$dmoment$dOmega,
                                               d2mu = out$d2moment$d2mu,
                                               d2Omega = out$d2moment$d2Omega,
                                               OmegaM1 = out$moment$OmegaM1.missing.pattern,
                                               missing.pattern = out$missing$pattern,
                                               unique.pattern = out$missing$unique.pattern,
                                               name.pattern = out$missing$name.pattern,
                                               grid.3varD1 = out$skeleton$grid.3varD1,
                                               grid.2meanD1.1varD1 = out$skeleton$grid.2meanD1.1varD1,
                                               grid.2meanD2.1meanD1 = out$skeleton$grid.2meanD2.1meanD1,
                                               grid.2varD2.1varD1 = out$skeleton$grid.2varD2.1varD1,
                                               name.param = out$skeleton$Uparam,
                                               leverage = out$leverage,
                                               n.cluster = out$cluster$n.cluster,
                                               weights = out$weights)

            ## delta method
            out$dVcov.param <- .dVcov.param(vcov.param = out$vcov.param,
                                            dInformation = out$dInformation,
                                            n.param = length(out$skeleton$Uparam),
                                            name.param = out$skeleton$Uparam)
        }

    ## *** dRvcov.param  (robust variance, analytic)
    if((derivative == "analytic") && dVcov.robust){

        out$dRvcov.param <- .dRvcov.param(score = out$score,
                                          hessian = out$hessian,
                                          vcov.param = out$vcov.param,
                                          dVcov.param = out$dVcov.param,
                                          n.param = length(out$skeleton$Uparam),
                                          name.param = out$skeleton$Uparam)
    }

    ## *** dVcov.param and dRvcov.param (numeric derivatives)
    if(derivative == "numeric" && (dVcov || dVcov.robust)){
        test.package <- try(requireNamespace("numDeriv"), silent = TRUE)
        if(inherits(test.package,"try-error")){
            stop("There is no package \'numDeriv\' \n",
                 "This package is necessary when argument \'numeric.derivative\' is TRUE \n")
        }

        ## range(out$score - .warper.numDev(value = out$param, object = object, type = "score"))
        ## range(out$hessian - .warper.numDev(value = out$param, object = object, type = "hessian"))
        ## range(out$vcov.param - .warper.numDev(value = out$param, object = object, type = "vcov.model"))

        param <- out$param
        name.param <- names(out$param)
        n.param <- length(param)
        n.cluster <- out$cluster$n.cluster
        object2 <- object
        object2$sCorrect <- out
        
        ## *** hessian
        num.hessian <- numDeriv::jacobian(.warper.numDev, x = param, object = object2, type = "score", method = "Richardson")

        out$hessian <- aperm(array(num.hessian, dim = c(n.cluster,n.param,n.param),
                                               dimnames = list(NULL, name.param, name.param)), perm = 3:1)

        ## ## *** dInformation
        num.information <- numDeriv::jacobian(.warper.numDev, x = param, object = object2, type = "information", method = "Richardson")

        out$dInformation <- array(num.information, dim = c(n.param,n.param,n.param),
                                              dimnames = list(name.param, name.param, name.param))

        ## *** dVcov.param
        num.dVcov.param <- numDeriv::jacobian(.warper.numDev, x = param, object = object2, type = "vcov.model", method = "Richardson")

        out$dVcov.param <- array(num.dVcov.param, dim = c(n.param,n.param,n.param),
                                             dimnames = list(name.param, name.param, name.param))


        ## *** dRvcov.param
        num.dRvcov.param <- numDeriv::jacobian(.warper.numDev, x = param, object = object2, type = "vcov.robust", method = "Richardson")

        out$dRvcov.param <- array(num.dRvcov.param, dim = c(n.param,n.param,n.param),
                                              dimnames = list(name.param, name.param, name.param))
    }        
    
        
    }

    ## ** export
    return(out)
}

## * moments2.lvmfit
#' @rdname moments2
#' @export
moments2.lvmfit <- moments2.lvm

## * .wraper.numDev (helper)
.warper.numDev <- function(value, object, type){ # x <- p.obj

    type <- match.arg(type, c("score","hessian","information","vcov.model","vcov.robust"))
            
    ## update moments and their derivatives (and also residuals! This makes a difference when computing the hessian)
    cM <- moments2(object, param = value,
                   initialize = FALSE, first.order = TRUE, second.order = (type == "hessian"),
                   usefit = TRUE, residuals = TRUE, leverage = FALSE)

    ## compute information matrix
    if(type %in% c("score","vcov.robust")){
        ss <- .score2(dmu = cM$dmoment$dmu,
                      dOmega = cM$dmoment$dOmega,
                      epsilon = cM$residuals,
                      OmegaM1 = cM$moment$OmegaM1.missing.pattern,
                      missing.pattern = object$sCorrect$missing$pattern,
                      unique.pattern = object$sCorrect$missing$unique.pattern,
                      name.pattern = object$sCorrect$missing$name.pattern,
                      name.param = object$sCorrect$skeleton$Uparam,
                      name.meanparam = object$sCorrect$skeleton$Uparam.mean,
                      name.varparam = object$sCorrect$skeleton$Uparam.var,
                      n.cluster = object$sCorrect$cluster$n.cluster,
                      weights = object$sCorrect$weights)
    }
    if(type=="hessian"){
        hh <- .hessian2(dmu = cM$dmoment$dmu,
                        dOmega = cM$dmoment$dOmega,
                        d2mu = cM$d2moment$d2mu,
                        d2Omega = cM$d2moment$d2Omega,
                        epsilon = cM$residuals,
                        OmegaM1 = cM$moment$OmegaM1.missing.pattern,
                        missing.pattern = object$sCorrect$missing$pattern,
                        unique.pattern = object$sCorrect$missing$unique.pattern,
                        name.pattern = object$sCorrect$missing$name.pattern,
                        grid.mean = object$sCorrect$skeleton$grid.dmoment$mean,
                        grid.var = object$sCorrect$skeleton$grid.dmoment$var,
                        grid.hybrid = object$sCorrect$skeleton$grid.dmoment$hybrid,
                        name.param = object$sCorrect$skeleton$Uparam,
                        leverage = object$sCorrect$leverage,
                        n.cluster = object$sCorrect$cluster$n.cluster,
                        weights = object$sCorrect$weights)
    }
    if(type %in% c("information","vcov.model","vcov.robust")){
        ## print(object$sCorrect$cluster$n.cluster - colSums(object$sCorrect$leverage))
        ii <- .information2(dmu = cM$dmoment$dmu,
                            dOmega = cM$dmoment$dOmega,
                            OmegaM1 = cM$moment$OmegaM1.missing.pattern,
                            missing.pattern = object$sCorrect$missing$pattern,
                            unique.pattern = object$sCorrect$missing$unique.pattern,
                            name.pattern = object$sCorrect$missing$name.pattern,
                            grid.mean = object$sCorrect$skeleton$grid.dmoment$mean,
                            grid.var = object$sCorrect$skeleton$grid.dmoment$var,
                            name.param = object$sCorrect$skeleton$Uparam,
                            leverage = object$sCorrect$leverage,
                            n.cluster = object$sCorrect$cluster$n.cluster,
                            weights = object$sCorrect$weights)
    }

    return(switch(type,
                  "score" = ss,
                  "hessian" = hh,
                  "information" = ii,
                  "vcov.model" = solve(ii),
                  "vcov.robust" = crossprod(ss %*% solve(ii)))
           )
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

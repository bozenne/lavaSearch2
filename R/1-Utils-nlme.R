### coef2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 15 2017 (17:29) 
## Version: 
## Last-Updated: jan 10 2018 (13:26) 
##           By: Brice Ozenne
##     Update #: 123
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .coef2
`.coef2` <-
    function(object) UseMethod(".coef2")

## * .coef2.gls
.coef2.gls <- function(object){

     ## *** mean parameters
    mean.coef <- coef(object)

    ## *** variance parameters
    if(!is.null(object$modelStruct$varStruct)){
        var.coef <- c(sigma2 = sigma(object),coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))^2
    }else{
        var.coef <- c(sigma2 = sigma(object)^2)
    }

    ## *** covariance parameters
    if(!is.null(object$modelStruct$corStruct)){
        cor.coef <- coef(object$modelStruct$corStruct, unconstrained = FALSE)
        names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
    }else{
        cor.coef <- NULL
    }

    p <- c(mean.coef, cor.coef, var.coef)
    attr(p, "mean.coef") <- names(mean.coef)
    attr(p, "var.coef") <- names(var.coef)
    attr(p, "cor.coef") <- names(cor.coef)
    return(p)
}




## * .coef2.lme
.coef2.lme <- function(object){

     ## *** mean parameters
    mean.coef <- fixef(object)

    ## *** variance parameters
    if(!is.null(object$modelStruct$varStruct)){
        var.coef <- c(sigma2 = sigma(object),coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))^2
    }else{
        var.coef <- c(sigma2 = sigma(object)^2)
    }

    ## *** random effect parameters
    random.coef <- as.double(getVarCov(object))    
    names(random.coef) <- paste0("ranCoef",1:length(random.coef))

     ## *** correlation parameters
    if(!is.null(object$modelStruct$corStruct)){
        cor.coef <- coef(object$modelStruct$corStruct, unconstrained = FALSE)
        names(cor.coef) <- paste0("corCoef",1:length(cor.coef))
    }else{
        cor.coef <- NULL
    }
    
    p <- c(mean.coef, cor.coef, var.coef, random.coef)
    attr(p, "mean.coef") <- names(mean.coef)
    attr(p, "var.coef") <- names(var.coef)
    attr(p, "cor.coef") <- names(cor.coef)
    attr(p, "ran.coef") <- names(random.coef)
    return(p)
}

## * .getFormula2
`.getFormula2` <-
    function(object) UseMethod(".getFormula2")

## * .getFormula2.gls
.getFormula2.gls <- function(object){
    return(evalInParentEnv(object$call$model))
}

## * .getFormula2.lme
.getFormula2.lme <- function(object){
    return(evalInParentEnv(object$call$fixed))
}


## * .getGroups2
`.getGroups2` <-
    function(object, ...) UseMethod(".getGroups2")

## * .getGroups2.gls
.getGroups2.gls <- function(object, cluster, data, ...){

    ### ** get cluster
    cluster2 <- as.numeric(nlme::getGroups(object))
    if(length(cluster2)==0){ ## no correlation
        if(length(cluster) == 1 && is.character(cluster)){
            cluster2 <- as.numeric(as.factor(data[[cluster]]))
        }else{
            if(length(cluster)!=NROW(data)){
                stop("length of cluster and data do not match \n")
            }
            cluster2 <- as.numeric(as.factor(cluster))
        }
    }
    n.cluster <- length(unique(cluster2))

    ### ** get outcome
    if(!is.null(object$modelStruct$varStruct)){
        name.endogenous <- attr(object$modelStruct$varStruct,"groupName")
        vec.rep <- attr(object$modelStruct$varStruct,"groups")        
    }else if(!is.null(object$modelStruct$corStruct)){        
        vec.rep0 <- unlist(attr(object$modelStruct$corStruct,"covariate"))
        if(is.factor(vec.rep0)){
            name.endogenous <- levels(vec.rep0)
        }else{
            name.endogenous <- unique(vec.rep0)
        }
        vec.rep <- vec.rep0[order(order(cluster2))]
    }else{
        vec.rep <- as.data.table(cluster2)[, rep := 1:.N, by = cluster2][["rep"]]
        # vec.rep <- rep("1",NROW(data))
        name.endogenous <- unique(vec.rep)
    }
    vec.rep <- as.numeric(factor(vec.rep, levels = name.endogenous))
    n.endogenous <- length(name.endogenous)
    
    ### ** convert observations from the vector format to the matrix format
    index.obs <- cluster2+(vec.rep-1)*n.cluster

    ### ** export
    name.endogenous <- as.character(name.endogenous)

    return(list(cluster = cluster2,
                n.cluster = n.cluster,
                endogenous = vec.rep,
                name.endogenous = name.endogenous,
                n.endogenous = n.endogenous,
                index.obs = index.obs))
}

## * .getGroups2.lme
.getGroups2.lme <- function(object, ...){

    ## ** get cluster
    if(NCOL(object$groups)!=1){
        stop("cannot only handle one random effect \n")
    }
    cluster <- as.numeric(nlme::getGroups(object))
    n.cluster <- length(unique(cluster))

    ## ** get outcome
    if(!is.null(object$modelStruct$varStruct)){
        name.endogenous <- attr(object$modelStruct$varStruct,"groupName")
        vec.rep <- attr(object$modelStruct$varStruct,"groups")        
    }else{
        if(!is.null(object$modelStruct$corStruct)){        
            vec.rep0 <- unlist(attr(object$modelStruct$corStruct,"covariate"))
        }else{
            vec.rep0 <- unlist(tapply(cluster,cluster,function(x){1:length(x)}))
        }
    
        if(is.factor(vec.rep0)){
            name.endogenous <- levels(vec.rep0)
        }else{
            name.endogenous <- unique(vec.rep0)
        }
        vec.rep <- vec.rep0[order(order(cluster))]
    }
    vec.rep <- as.numeric(factor(vec.rep, levels = name.endogenous))
    n.endogenous <- length(name.endogenous)
    
    ## ** convert observations from the vector format to the matrix format
    index.obs <- cluster+(vec.rep-1)*n.cluster

### ** export
    name.endogenous <- as.character(name.endogenous)
    
    return(list(cluster = cluster,
                n.cluster = n.cluster,
                endogenous = vec.rep,
                name.endogenous = name.endogenous,
                n.endogenous = n.endogenous,
                index.obs = index.obs))
}

## * .getVarCov2
#' @title Reconstruct the marginal variance covariance matrix from a nlme model
#' @description Reconstruct the marginal variance covariance matrix from a nlme model
#' @name getVarCov2
#'
#' @param object a gls or lme object
#' @param ... arguments to be passed to lower level methods.
#'  
#' @details The marginal variance covariance matrix for gls model is of the form:
#' \deqn{
#' \Sigma = \left[ \begin{array}{ccc}
#' \sigma^2 & \sigma^2 \sigma_2 \rho_{1,2} & \sigma^2 \sigma_3 \rho_{1,3} \\
#' & \sigma^2 \sigma_2^2 & \sigma^2 \sigma_3 \rho_{1,3} \\
#' & & \sigma^2 \sigma_3^2 \\
#' \end{array} \right]
#' }
#'
#' The marginal variance covariance matrix for lme model is of the form:
#' 
`.getVarCov2` <-
    function(object, ...) UseMethod(".getVarCov2")

## * .getVarCov2.gls
.getVarCov2.gls <- function(object, param, attr.param,
                            endogenous, name.endogenous, n.endogenous,
                            cluster, n.cluster){

    var.coef <- param[attr.param$var.coef]
    cor.coef <- param[attr.param$cor.coef]

    ## ** Diagonal terms
    name.other <- setdiff(names(var.coef),"sigma2")
    if(length(name.other)>0){            
        sigma2.base <- setNames(var.coef["sigma2"]*c(1,var.coef[name.other]), name.endogenous)            
    }else{
        sigma2.base <- setNames(rep(var.coef["sigma2"],n.endogenous), name.endogenous)
    }
    template <- diag(as.double(sigma2.base),
                     nrow = n.endogenous, ncol = n.endogenous)
    
    ## ** Extra-diagonal terms
    if(length(cor.coef)>0){
        index.lower <- which(lower.tri(template))
        index.lower.arr <- which(lower.tri(template),arr.ind = TRUE)
        vec.sigma.tempo <- apply(index.lower.arr,1,function(x){prod(sqrt(sigma2.base[x]))})        
        template[index.lower] <- cor.coef*vec.sigma.tempo
        template <- symmetrize(template)
    }    
    dimnames(template) <- list(name.endogenous, name.endogenous)
    
    ## ** Index for each cluster    
    ls.indexOmega <- tapply(endogenous, cluster, function(iRep){iRep})

    ## ** export
    return(list(Omega = template,
                ls.indexOmega = ls.indexOmega))
}

## * .getVarCov2.lme
.getVarCov2.lme <- function(object, param, attr.param,
                            endogenous, name.endogenous, n.endogenous,
                            cluster, n.cluster){

    ## ** prepare with gls
    out <- .getVarCov2.gls(object, param = param, attr.param = attr.param,
                    endogenous = endogenous, name.endogenous = name.endogenous, n.endogenous = n.endogenous,
                    cluster = cluster, n.cluster = n.cluster)

    ## ** add contribution of the random effect
    ran.coef <- param[attr.param$ran.coef]
    out$Omega <- out$Omega + ran.coef

    ## ** export
    return(out)
    
}




##----------------------------------------------------------------------
### coef2.R ends here

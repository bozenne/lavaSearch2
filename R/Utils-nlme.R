### coef2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 15 2017 (17:29) 
## Version: 
## Last-Updated: nov 15 2017 (18:40) 
##           By: Brice Ozenne
##     Update #: 51
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

## * .coef.lme
.coef2.lme <- function(object){

     ## *** mean parameters
    mean.coef <- fixef(object)

    ## *** variance parameters
    if(!is.null(object$modelStruct$varStruct)){
        var.coef <- c(sigma2 = sigma(object),coef(object$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))^2
    }else{
        var.coef <- c(sigma2 = sigma(object)^2)
    }

    ## *** covariance parameters
    cor.coef <- as.double(getVarCov(object))    
    names(cor.coef) <- paste0("corCoef",1:length(cor.coef))

    p <- c(mean.coef, cor.coef, var.coef)
    attr(p, "mean.coef") <- names(mean.coef)
    attr(p, "var.coef") <- names(var.coef)
    attr(p, "cor.coef") <- names(cor.coef)
    return(p)
}

## * .coef.gls
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



##----------------------------------------------------------------------
### coef2.R ends here


## * .getVarCov2
`.getVarCov2` <-
    function(object, ...) UseMethod(".getVarCov2")

## * .getVarCov2.lme
.getVarCov2.lme <- function(object, var.coef, cor.coef){

    
}

## * .getVarCov2.gls
.getVarCov2.gls <- function(object, var.coef, cor.coef,
                            n.coef,
                            endogenous, name.endogenous, n.endogenous,
                            cluster, n.cluster){

    ## ** Diagonal terms
    name.other <- setdiff(names(var.coef),"sigma2")
    if(length(name.other)){            
        sigma2.base <- setNames(var.coef["sigma2"]*c(1,var.coef[name.other]), name.endogenous)            
    }else{
        sigma2.base <- setNames(rep(var.coef["sigma2"],n.endogenous), name.endogenous)
    }
    template <- diag(sigma2.base)
    
    ## ** Extra-diagonal terms
    if(length(cor.coef)>0){
        index.lower <- which(lower.tri(template))
        index.lower.arr <- which(lower.tri(template),arr.ind = TRUE)
        vec.sigma.tempo <- apply(index.lower.arr,1,function(x){prod(sqrt(sigma2.base[x]))})        
        template[index.lower] <- cor.coef*vec.sigma.tempo
    }

    
    ## ** Individual matrix
    template <- symmetrize(template)
    dimnames(template) <- list(name.endogenous, name.endogenous)
    
    Omega <- tapply(endogenous, cluster, function(iRep){
        return(list(template[iRep,iRep,drop=FALSE]))
    })

    return(Omega)
}

## * .getGroups2
`.getGroups2` <-
    function(object, ...) UseMethod(".getGroups2")

## * .getGroups2.lme
.getGroups2.lme <- function(object, ...){

    ## ** get cluster
    cluster2 <- as.numeric(nlme::getGroups(object))

    ### ** get outcome
    return(cluster2)
}

## * .getGroups2.gls
.getGroups2.gls <- function(object, cluster, data, ...){

    ### ** get cluster
    cluster2 <- as.numeric(nlme::getGroups(object))
    if(is.null(cluster2)){ ## no correlation
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

    ## ** get outcome
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
        vec.rep <- rep("1",n.cluster2)
        name.endogenous <- "1"
    }
    n.endogenous <- length(name.endogenous)
    
    ## ** convert observations from the vector format to the matrix format
    index.obs <- cluster2+(match(vec.rep,name.endogenous)-1)*n.cluster

    
    return(list(cluster = cluster2,
                n.cluster = n.cluster,
                endogenous = vec.rep,
                name.endogenous = name.endogenous,
                n.endogenous = n.endogenous,
                index.obs = index.obs))
}


##----------------------------------------------------------------------
### coef2.R ends here

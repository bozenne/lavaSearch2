### coefType.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (14:38) 
## Version: 
## last-updated: okt 13 2017 (16:42) 
##           By: Brice Ozenne
##     Update #: 77
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - coefType
#' @title Extract the specific coefficient names or positions in a LVM
#' @description Extract the specific coefficient names or positions in a LVM
#' 
#' @name coefType
#' 
#' @param x a lvm model or a fitted lvm model 
#' @param value should the name of the coefficient be returned? Else return the coefficients
#' @param detailed should the type of parameter be returned as a mathematical symbol? Otherwise plain english.
#' @param keep.var should the variance parameters be output?
#' @param level level argument of \code{lava::coef}
#' @param ... arguments to be passed to \code{lava::coef}
#'
#' @details A lvm can be written as a measurement model:
#' \deqn{Y_i = \nu + \Lambda \eta_i + K X_i + \epsilon_i}
#' and a structural model:
#' \deqn{\eta_i = \alpha + B \eta_i + \Gamma X_i + \zeta_i}
#' where \eqn{\Psi}   is the variance covariance matrix of the residuals \eqn{\zeta} \cr
#' and   \eqn{\Sigma} is the variance covariance matrix of the residuals \eqn{\epsilon}. \cr \cr
#'
#' \code{coefType} either returns the latin/greek letter corresponding to the coefficients
#' or it groups them:
#' \itemize{
#' \item intercept: \eqn{\nu} and \eqn{\alpha}.
#' \item regression: \eqn{\Lambda}, \eqn{K}, \eqn{B}, and \eqn{\Gamma}.
#' \item covariance: extra-diagonal terms of \eqn{\Sigma} and \eqn{\Psi}.
#' \item variance: diagonal of \eqn{\Sigma} and \eqn{\Psi}.
#' }
#' 
#' @examples 
#' #### regression ####
#' m <- lvm(Y~X1+X2)
#' e <- estimate(m, sim(m, 1e2))
#' 
#' coefType(m)
#' coefType(e)
#' coefType(e, level = -1)
#'
#' coefCov(m)
#' coefCov(m, value = TRUE)
#'
#' coefCov(m, keep.var = TRUE)
#' coefCov(m, value = TRUE, keep.var = TRUE)
#'
#' coefIndexModel(m)
#' coefIndexModel(e)
#' 
#' coefIntercept(m)
#' coefIntercept(m, value = TRUE)
#'
#' coefReg(m)
#' coefReg(m, value = TRUE)
#' 
#' #### LVM ####
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x1+x2
#' latent(m) <- ~u
#' covariance(m) <- y1~y2
#'
#' m.Sim <- m
#' categorical(m.Sim, labels = c("a","b","c")) <- ~x2
#' e <- estimate(m, sim(m.Sim, 1e2))
#' 
#' coefType(m)
#' coefType(e)
#' coefType(e, level = -1)
#'
#' coefCov(m)
#' coefCov(m, value = TRUE)#' 
#'
#' coefCov(m, keep.var = TRUE)
#' coefCov(m, value = TRUE, keep.var = TRUE)
#' 
#' coefExtra(m)
#'
#' coefIndexModel(m)
#' coefIndexModel(e)
#' 
#' categorical(m, K = 3) <- "X1"
#' coefExtra(m)
#' coefExtra(m, value = TRUE)
#'
#' coefIntercept(m)
#' coefIntercept(m, value = TRUE)
#' coefIntercept(e)
#'
#' coefReg(e, value = TRUE)
#' coefReg(e, level = -1, value = TRUE)
#'
#' #### multigroup #####' 
#' eG <- estimate(list(m,m), list(sim(m, 1e2),sim(m, 1e2)))
#' coefType(eG)
#'
#' @export
`coefType` <-
  function(x,...) UseMethod("coefType")


## * method - coefType
## ** method - coefType.lvm
#' @rdname coefType
#' @export
coefType.lvm <- function(x, detailed = FALSE, ...){ 

### ** extract all coef
    names.coef <- coef(x, ...)
    index.coef <- lava::index(x)

### ** find coef class
    type <- setNames(character(length = length(names.coef)), names.coef)
    type[index.coef$parBelongsTo$mean] <- "intercept"
    type[index.coef$parBelongsTo$reg] <- "regression"
    type[index.coef$parBelongsTo$cov] <- "covariance"
    type[index.coef$parBelongsTo$epar] <- "extra"
  
    type[names(names.coef) %in% diag(APmatrix(x)$P)] <- "variance"

### ** rename as
    if(detailed){
        ls.links <- initVarLinks(names.coef)
        index.loading <- setdiff(which(ls.links$var2 %in% latent(x)),
                                 which(type %in% c("covariance","variance")))
        if(length(index.loading)>0){
            type[index.loading] <- "loading"
        }
        
        index.measurement <- which(ls.links$var1 %in% endogenous(x))
        if(length(index.measurement)>0){
            type[index.measurement] <- as.character(factor(type[index.measurement],
                                                           levels = c("intercept","regression","loading","covariance","variance"),
                                                           labels = c("nu","K","Lambda","Sigma_cov","Sigma_var")))
        }
        index.structural <- setdiff(1:length(type),index.measurement)
        if(length(index.structural)>0){
            type[index.structural] <- as.character(factor(type[index.structural],
                                                          levels = c("intercept","regression","loading","covariance","variance"),
                                                          labels = c("alpha","Gamma","B","Psi_cov","Psi_var")))
        }
        
    }
### ** export
    return(type)
}

## ** method - coefType.lvmfit
#' @rdname coefType
#' @export
coefType.lvmfit <- function(x, level = 9, index.model = FALSE, ...){ 
  
### ** extract all coefficients
    extractCoef <- coef(x, level = level, ...)
    if(is.matrix(extractCoef)){
        names.allCoef <- rownames(extractCoef)
    }else{
        names.allCoef <- names(extractCoef)
    }
    
    res <- setNames(rep(NA, length = length(names.allCoef)), names.allCoef)
    attribute <- setNames(rep(TRUE, length = length(names.allCoef)), names.allCoef)

### ** find type of the coefficients in the original model
    type <- coefType(x$model0, ...)        
    res[names(type)] <- type
    attribute[names(type)] <- FALSE

### ** find the missing coefficients which are in fact fixed as reference
    originalX <- evalInParentEnv(x$call$x)
    if("lvm" %in% class(originalX)){
        reftype <- coefType(originalX, ...)
        res[attribute==TRUE] <- reftype[names(attribute)[attribute==TRUE]]
    }   

    attr(res,"reference") <- attribute
    if(index.model){
        attr(type, "index.model") <- rep(1, length(type))
    }
    
    #### export
    return(res)
}
# }}}

## ** method - coefType.multigroup
#' @rdname coefType
#' @export
coefType.multigroup <- function(x, ...){ 
  n.model <- length(x$lvm)
  
  ## new coef names
  allCoef <- x$name
  n.allCoef <- length(allCoef)
  index.AllCoef <- x$coef
  
  type_tempo <- NULL
  type <- setNames(rep("", n.allCoef), allCoef)
  
  ## old coef names
  indexCoef.old <- x$coef.idx
  for(iModel in 1:n.model){ # iModel <- 1
    
      if(!is.null(x$meanposN[[iModel]])){
          indexCoef.old[[iModel]] <- c(1,indexCoef.old[[iModel]]+1)
      }

      type_tempo <- coefType(x$lvm[[iModel]],...)      
      type[index.AllCoef[[iModel]]] <- type_tempo[indexCoef.old[[iModel]]]
  }
  
  #### export
  return(type)
}

## * method - coefCov
#' @rdname coefType
#' @export
`coefCov` <-
  function(x,...) UseMethod("coefCov")

#' @rdname coefType
#' @export
coefCov.lvm <- function(x, value = FALSE, keep.var = FALSE, ...){

    res <- retainType(type = coefType(x, ...),
                      validType = c("covariance", if(keep.var){"variance"}else{NULL}),
                      value = value)

    return(res)
}

#' @rdname coefType
#' @export
coefCov.lvmfit <- coefCov.lvm

#' @rdname coefType
#' @export
coefCov.multigroup <- coefCov.lvm

## * method - coefExtra

#' @rdname coefType
#' @export
`coefExtra` <-
  function(x,...) UseMethod("coefExtra")

#' @rdname coefType
#' @export
coefExtra.lvm <- function(x, value = FALSE, ...){ 

    res <- retainType(type = coefType(x, ...),
                      validType = "extra",
                      value = value) 
    
    return(res)    
}

#' @rdname coefType
#' @export
coefExtra.lvmfit <- coefExtra.lvm

#' @rdname coefType
#' @export
coefExtra.multigroup <- coefExtra.lvm

## * method - coefIndexModel
#' @rdname coefType
#' @export
`coefIndexModel` <-
  function(x,...) UseMethod("coefIndexModel")

#' @rdname coefType
#' @export
coefIndexModel.lvm <- function(x, ...){
    name.coef <- coef(x)
    index <- rep(1, length(name.coef))
    names(index) <- name.coef
    return(index)
}

#' @rdname coefType
#' @export
coefIndexModel.lvmfit <- function(x, ...){
    name.coef <- names(coef(x))
    index <- rep(1, length(name.coef))
    names(index) <- name.coef
    return(index)
}

#' @rdname coefType
#' @export
coefIndexModel.multigroup <- function(x, ...){
    n.model <- length(x$lvm)

    ## new coef names
    allCoef <- x$name
    n.allCoef <- length(allCoef)
    index.AllCoef <- x$coef
  
    index <- setNames(rep(NA, n.allCoef), allCoef)
  
    for(iModel in 1:n.model){ # iModel <- 1
        index[index.AllCoef[[iModel]]] <- iModel
    }
  
    #### export
    return(index)
}
  
# }}}

## * method - coefIntercept

#' @rdname coefType
#' @export
`coefIntercept` <-
  function(x,...) UseMethod("coefIntercept")

#' @rdname coefType
#' @export
coefIntercept.lvm <- function(x, value = FALSE, ...){ 

    res <- retainType(type = coefType(x, ...),
                      validType = "intercept",
                      value = value)

    return(res)
}

#' @rdname coefType
#' @export
coefIntercept.lvmfit <- coefIntercept.lvm

#' @rdname coefType
#' @export
coefIntercept.multigroup <- coefIntercept.lvm

## * method - coefRef
#' @rdname coefType
#' @export
`coefRef` <-
  function(x,...) UseMethod("coefRef")

#' @rdname coefType
#' @export
coefRef.lvmfit <- function(x, value = FALSE, ...){
    
    res <- retainType(type = attr(coefType(x, ...), "reference"),
                      validType = TRUE,
                      value = value)

    return(res)    
}

# }}}

#' @rdname coefType
#' @export
`coefReg` <-
  function(x,...) UseMethod("coefReg")

#' @rdname coefType
#' @export
coefReg.lvm <- function(x, value = FALSE, ...){
    
     res <- retainType(type = coefType(x, ...),
                      validType = "regression",
                      value = value)

     return(res)
}

#' @rdname coefType
#' @export
coefReg.lvmfit <- coefReg.lvm

#' @rdname coefType
#' @export
coefReg.multigroup <- coefReg.lvm

## * method - coefVar

#' @rdname coefType
#' @export
`coefVar` <-
  function(x,...) UseMethod("coefVar")

#' @rdname coefType
#' @export
coefVar.lvm <- function(x, value = FALSE, ...){ 

    res <- retainType(type = coefType(x, ...),
                      validType = "variance",
                      value = value)

    return(res)
}

#' @rdname coefType
#' @export
coefVar.lvmfit <- coefVar.lvm

#' @rdname coefType
#' @export
coefVar.multigroup <- coefVar.lvm


#----------------------------------------------------------------------
### coefType.R ends here

## * retainType  (need for coefCov/Latent/Ref)
retainType <- function(type, validType, value){
  index.var <- which(type %in% validType)
  
  if(length(index.var)>0){
      if(value){
          return(names(type)[index.var])
      }else{
          return(index.var)
      }
  }else{
      return(NULL)
  }
}

## * APmatrix (need for coefType)
APmatrix <- function(x){ # borrowed from coef.lvmfit
 
  names2.coef <- names(coef(x))
  if (is.null(x$control$meanstructure)){
    meanstructure <- TRUE
  } else {
    meanstructure <- x$control$meanstructure
  }
  npar <- lava::index(x)$npar
  npar.mean <- lava::index(x)$npar.mean*meanstructure
  npar.ex <- lava::index(x)$npar.ex
  
  if (inherits(x,"lvm.missing")) {
    if (length(x$cc)==0) {## No complete cases
      coefs <- coef(x$estimate)
      c1 <- coef(Model(x),mean=TRUE,fix=FALSE)
      c1. <- coef(Model(x),mean=FALSE,fix=FALSE)
      myorder <- match(c1,names(coefs))
      myorder.reg <- order(na.omit(match(names(coefs),c1.)))
      myorder.extra <- c()
      ##mp <-effect modelPar(x,seq_len(npar+npar.mean+npar.ex))
      ## mp <- modelPar(x,seq_len(npar+npar.mean+npar.ex))
      ## myorder <- c(mp$meanpar,mp$p)
      ## myorder.reg <- seq_len(length(mp$p))
      ## myorder.extra <- mp$p2
    } else {
      myorder <- na.omit(modelPar(x$multigroup,seq_len(npar+npar.mean))$p[[x$cc]])
      myorder.reg <- na.omit(modelPar(x$multigroup,seq_len(npar))$p[[x$cc]])
      myorder.extra <- seq_len(lava::index(x)$npar.ex)+length(myorder)
      myorder <- c(myorder,myorder.extra)
    }
  } else {
    myorder <- seq_len(npar+npar.mean)
    myorder.reg <- seq_len(npar)
    myorder.extra <- seq_len(lava::index(x)$npar.ex)+length(myorder)
    myorder <- c(myorder,myorder.extra)
  }
  ## myorder <- seq_len(npar+npar.mean)
  ## myorder.reg <- seq_len(npar)
  ## myorder.extra <- seq_len(lava::index(x)$npar.ex)+length(myorder)
  ## myorder <- c(myorder,myorder.extra)
  
  myparnames <- paste0("p",seq_len(npar+npar.ex))[myorder.reg]
  return(lava_matrices.lvm(Model(x), myparnames))
  
}

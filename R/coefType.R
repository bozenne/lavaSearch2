### coefType.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (14:38) 
## Version: 
## last-updated: nov 28 2019 (11:35) 
##           By: Brice Ozenne
##     Update #: 718
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - coefType
#' @title Extract the Type of Each Coefficient
#' @description Extract the type of each coefficient of a \code{lvm} object.
#' @name coefType
#' 
#' @param object a \code{lvm} or \code{lvmfit} object. 
#' @param data [data.frame, optional] the dataset. Help to identify the categorical variables.
#' @param as.lava [logical] export the type of coefficients mimicking \code{lava:::coef}.
#' @param attr.type [logical] add as an attribute the list containing parameters by type.
#' @param ... arguments to be passed to \code{lava::coef}
#'
#' @details A lvm can be written as a measurement model:
#' \deqn{Y_i = \nu + \Lambda \eta_i + K X_i + \epsilon_i}
#' and a structural model:
#' \deqn{\eta_i = \alpha + B \eta_i + \Gamma X_i + \zeta_i}
#' where \eqn{\Psi}   is the variance covariance matrix of the residuals \eqn{\zeta} \cr
#' and   \eqn{\Sigma} is the variance covariance matrix of the residuals \eqn{\epsilon}. \cr \cr
#'
#' \code{coefType} either returns the Latin/Greek letter corresponding to the coefficients
#' or it groups them:
#' \itemize{
#' \item intercept: \eqn{\nu} and \eqn{\alpha}.
#' \item regression: \eqn{\Lambda}, \eqn{K}, \eqn{B}, and \eqn{\Gamma}.
#' \item covariance: extra-diagonal terms of \eqn{\Sigma} and \eqn{\Psi}.
#' \item variance: diagonal of \eqn{\Sigma} and \eqn{\Psi}.
#' }
#'
#' A link denotes a relationship between two variables.
#' The coefficient are used to represent the strength of the association between two variable, i.e. the strength of a link.
#' A coefficient may corresponds to the strength of one or several link.
#'
#' @return \code{coefType} returns a \code{data.frame} when \code{as.lava=FALSE}:
#' \itemize{
#' \item name: name of the link
#' \item Y: outcome variable
#' \item X: regression variable in the design matrix (could be a transformation of the original variables, e.g. dichotomization).
#' \item data: original variable
#' \item type: type of link
#' \item value: if TRUE, the value of the link is set and not estimated.
#' \item marginal: if TRUE, the value of the link does not impact the estimation.
#' \item detail: a more detailed description of the type of link (see the details section)
#' \item lava: name of the coefficient in lava
#' }
#' When \code{as.lava=TRUE}, \code{coefType} returns a named vector containing the type of each coefficient.
#' 
#' @examples 
#' #### regression ####
#' m <- lvm(Y~X1+X2)
#' e <- estimate(m, lava::sim(m, 1e2))
#' 
#' coefType(m)
#' coefType(e)
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
#' e <- estimate(m, lava::sim(m.Sim, 1e2))
#' 
#' coefType(m)
#' coefType(e)
#'
#' ## additional categorical variables 
#' categorical(m, labels = as.character(1:3)) <- "X1"
#'
#' coefType(m, as.lava = FALSE)
#'
#' #### LVM with constrains ####
#' m <- lvm(c(Y1~0+1*eta1,Y2~0+1*eta1,Y3~0+1*eta1,
#'           Z1~0+1*eta2,Z2~0+1*eta2,Z3~0+1*eta2))
#' latent(m) <- ~eta1 + eta2
#' e <- estimate(m, lava::sim(m,1e2))
#' 
#' coefType(m)
#' coefType(e)
#' 
#' #### multigroup ####
#' m <- lvm(Y~X1+X2)
#' eG <- estimate(list(m,m), list(lava::sim(m, 1e2), lava::sim(m, 1e2)))
#' coefType(eG)
#'
#' @concept extractor
#' @export
`coefType` <-
  function(object, as.lava, ...) UseMethod("coefType")

## * coefType.lm
#' @rdname coefType
#' @export
coefType.lm <- function(object, data = NULL, indexOmega = NULL, coef2 = NULL, attr.type = FALSE, ...){

    ## ** extract all coef
    if(is.null(data)){
        data <- extractData(object)
    }
    if(is.null(coef2)){
        coef2 <- coef2(object)
    }
    if(is.null(indexOmega)){
        indexOmega <- .getIndexOmega(object, data = data)
    }

    formula.object <- formula(object)
    object.terms <- terms(formula.object)
    name.coef <- names(coef2)
    name.endogenous <- endogenous(object, format = "long")
    name.endogenousW <- endogenous(object, format = "wide")
    n.endogenous <- length(name.endogenousW)
    name.exogenous <- exogenous(object)
    out <- NULL

    ## ** intercept
    n.intercept <- attr(object.terms, "intercept")
    if(n.intercept>0){
        out <- rbind(out,
                     data.frame(name = name.endogenousW,
                                Y = name.endogenousW,
                                X = as.character(NA),
                                data = as.character(NA),
                                type = "intercept",
                                value = as.double(coef2["(Intercept)"]),
                                param = "(Intercept)",
                                marignal = FALSE,
                                factitious = FALSE,
                                level = as.character(NA),
                                externalLink = as.character(NA),
                                originalLink = name.endogenousW,
                                detail = "nu",
                                lava = as.character(NA),
                                stringsAsFactors = FALSE)
                     )
    }

    ## ** regression
    regcoef <- setdiff(attr(coef2, "mean.coef"), "(Intercept)")
    n.regcoef <- length(regcoef)
    if(n.regcoef>0){
        
        ## ## rebuild X in presence of categorical variables
        ## if(inherits(object,"lm")){
        ##     factor.var <- names(object$xlevels)
        ##     if(any(sapply(object$contrast,function(iC){iC != "contr.treatment"}))){
        ##         stop("coefType only implemented for contr.treatment contrasts \n")
        ##     }
        ##     xlevels <- lapply(object$xlevels, function(iC){iC[-1]})
        ## }else{
        ##     factor.var <- names(object$contrasts)
        ##     xlevels <- lapply(object$contrasts, function(iC){ ## iC <- object$contrasts[[1]]
        ##         if(any(colSums(abs(iC))!=1)){
        ##             stop("coefType only implement for contr.treatment contrasts \n")
        ##         }
        ##         return(rownames(iC)[apply(iC==1,2,which)])
        ##     })
        ## }

        X <- model.matrix(formula.object,data)
        data.var <- setNames(attr(object.terms, "term.labels")[attr(X,"assign")[attr(X,"assign")!=0]],
                             setdiff(attr(coef2,"mean.coef"),"(Intercept)"))

        for(iE in 1:n.endogenous){
            index.endpoint <- which(indexOmega==iE)

            iTest.0 <- apply(X[index.endpoint,],2, function(iCol){any(iCol!=0)})
            iName.exogenous <- setdiff(names(iTest.0[iTest.0]),"(Intercept)")
            if(length(iName.exogenous)==0){
                next
            }
            iData.var <- data.var[iName.exogenous]
        
            out <- rbind(out,
                         data.frame(name = iName.exogenous,
                                    Y = name.endogenousW[iE],
                                    X = iName.exogenous,
                                    data = iData.var,
                                    type = "regression",
                                    value = as.double(coef2[iName.exogenous]),
                                    param = iName.exogenous,
                                    marignal = FALSE,
                                    factitious = FALSE,
                                    level = as.character(NA),
                                    externalLink = as.character(NA),
                                    originalLink = paste0(name.endogenousW[iE],"~",iName.exogenous),
                                    detail = "K",
                                    lava = as.character(NA),
                                    stringsAsFactors = FALSE)
                         )
        }
    }

    ## ** random effect
    if(!is.null(object$modelStruct$reStruct)){
        if(any("eta" %in% manifest(object))){
            stop("The name \"eta\" is used internally. Variables used in the model should not have this name.")
        }

        df.lambda <- data.frame(name = paste0(name.endogenousW,"~eta"),
                                Y = name.endogenousW,
                                X = "eta",
                                data = as.character(NA),
                                type = "variance",
                                value = 1,
                                param = as.character(NA),
                                marignal = FALSE,
                                factitious = FALSE,
                                level = as.character(NA),
                                externalLink = as.character(NA),
                                originalLink = paste0(name.endogenousW,"~eta"),
                                detail = "Lambda",
                                lava = as.character(NA),
                                stringsAsFactors = FALSE)
        df.eta <- data.frame(name = "eta~eta",
                             Y = "eta",
                             X = "eta",
                             data = as.character(NA),
                             type = "variance",
                             value = as.double(coef2[attr(coef2,"ran.coef")]),
                             param = attr(coef2,"ran.coef"),
                             marignal = FALSE,
                             factitious = FALSE,
                             level = as.character(NA),
                             externalLink = as.character(NA),
                             originalLink = "eta~eta",
                             detail = "Psi_var",
                             lava = as.character("r1"),
                             stringsAsFactors = FALSE)
        out <- rbind(out, df.lambda, df.eta)
    }

    ## ** correlation structure
    if(!is.null(object$modelStruct$corStruct)){
        cor.coef <- attr(coef2, "cor.coef")
        cor.combi <- .combination(name.endogenousW,name.endogenousW)
        cor.combi <- cor.combi[cor.combi[,1]!=cor.combi[,2],]
        
        out <- rbind(out,
                     data.frame(name = paste0(cor.combi[,1],"~~",cor.combi[,2]),
                                Y = cor.combi[,1],
                                X = cor.combi[,2],
                                data = name.endogenous,
                                type = "covariance",
                                value = as.double(coef2[cor.coef]),
                                param = cor.coef,
                                marignal = FALSE,
                                factitious = FALSE,
                                level = as.character(NA),
                                externalLink = as.character(NA),
                                originalLink = paste0(cor.combi[,1],"~~",cor.combi[,2]),
                                detail = "cor",
                                lava = as.character(NA),
                                stringsAsFactors = FALSE)
                     )
    }
    
    ## ** variance structure
    out <- rbind(out,
                 data.frame(name = paste0(name.endogenousW,"~~",name.endogenousW),
                            Y = name.endogenousW,
                            X = name.endogenousW,
                            data = name.endogenous,
                            type = "variance",
                            value = as.double(coef2["sigma2"]),
                            param = "sigma2",
                            marignal = FALSE,
                            factitious = FALSE,
                            level = as.character(NA),
                            externalLink = as.character(NA),
                            originalLink = paste0(name.endogenousW,"~~",name.endogenousW),
                            detail = if(is.null(object$modelStruct$varStruct) && is.null(object$modelStruct$corStruct)){"Sigma_var"}else{"sigma2"},
                            lava = as.character(NA),
                            stringsAsFactors = FALSE)
                 )

    if(!is.null(object$modelStruct$varStruct)){
        k.coef <- setdiff(attr(coef2, "var.coef"), "sigma2")

        out <- rbind(out,
                     data.frame(name = paste0(paste0(name.endogenous,k.coef),"~~",paste0(name.endogenous,k.coef)),
                                Y = paste0(name.endogenous,k.coef),
                                X = paste0(name.endogenous,k.coef),
                                data = name.endogenous,
                                type = "variance",
                                value = as.double(coef2[k.coef]),
                                param = k.coef,
                                marignal = FALSE,
                                factitious = FALSE,
                                level = as.character(NA),
                                externalLink = as.character(NA),
                                originalLink = paste0(paste0(name.endogenous,k.coef),"~~",paste0(name.endogenous,k.coef)),
                                detail = "sigma2k",
                                lava = as.character(NA),
                                stringsAsFactors = FALSE)
                     )
    }
    
    ## ** export
    test.param <- !duplicated(out$param)
    out$lava[test.param] <- paste0("p",cumsum(test.param[test.param]))
    rownames(out) <- NULL

    out$detail <- factor(out$detail,
                         levels = c("nu","alpha","K","Gamma","Lambda","B","Sigma_var","Sigma_cov","sigma2","sigma2k","cor","Psi_var","Psi_cov",NA))
    out <- out[order(out$detail,out$param),]

    if(attr.type){
        iTempo <- detailType(out, endogenous = endogenous(object), latent = latent(object))
        attr(out, "byType") <- iTempo$byType
        attr(out, "byTypebyEndo") <- iTempo$byTypebyEndo
    }
    return(out)
}

## * coefType.gls
#' @rdname coefType
#' @export
coefType.gls <- coefType.lm

## * coefType.lme
#' @rdname coefType
#' @export
coefType.lme <- coefType.lm

## * coefType.lvm
#' @rdname coefType
#' @export
coefType.lvm <- function(object, as.lava = TRUE, data = NULL, attr.type = FALSE, ...){

    externalLink <- type <- NULL ## [:for CRAN check] subset
    
    ## ** extract all coef
    index.all <- which(!is.na(object$M), arr.ind = FALSE)
    ls.name <- list()
    ls.X <- list()
    ls.Y <- list()
    ls.type <- list()
    ls.value <- list()
    ls.param <- list()
    ls.marginal <- list()

    ## ** intercept
    n.intercept <- length(object$mean)
    if(n.intercept>0){
        ls.name$intercept <- names(object$mean)    
    
        ls.Y$intercept <- ls.name$intercept
        ls.X$intercept <- rep(NA, n.intercept)    
        ls.type$intercept <- rep("intercept", n.intercept)
        ls.value$intercept <- lapply(object$mean, function(iP){if(is.numeric(iP)){iP}else{NA}})
        ls.param$intercept <- unlist(Map(function(iPar,iFix,iName){if(iFix){NA}else if(!is.na(iPar)){iPar} else {iName}},
                                         iPar = unlist(object$mean),
                                         iFix = !is.na(ls.value$intercept),
                                         iName = ls.name$intercept)
                                     )
        ls.marginal$intercept <-  ls.name$intercept %in% exogenous(object)
    }
    
    ## ** regression
    arrIndex.regression <- which(object$M==1, arr.ind = TRUE)
    index.regression <- which(object$M==1, arr.ind = FALSE)
    n.regression <- length(index.regression)
    if(n.regression>0){
    
        ls.Y$regression <- colnames(object$M)[arrIndex.regression[,"col"]]
        ls.X$regression <- rownames(object$M)[arrIndex.regression[,"row"]]
        ls.name$regression <- paste0(ls.Y$regression,
                                     lava.options()$symbols[1],
                                     ls.X$regression)
    
        ls.type$regression <- rep("regression", n.regression)    
        ls.value$regression <- object$fix[index.regression]
        ls.param$regression <- unlist(Map(function(iPar,iFix,iName){if(iFix){NA}else if(!is.na(iPar)){iPar} else {iName}},
                                          iPar = object$par[index.regression],
                                          iFix = !is.na(ls.value$regression),
                                          iName = ls.name$regression)
                                      )
        ls.marginal$regression <- rep(FALSE,n.regression)
    }

    ## ** covariance
    M.cov <- object$cov
    M.cov[upper.tri(M.cov)] <- 0
    
    arrIndex.vcov <- which(M.cov==1, arr.ind = TRUE)
    index.vcov <- which(M.cov==1, arr.ind = FALSE)
    n.vcov <- length(index.vcov)
    if(n.vcov>0){
        
        Y.vcov <- colnames(object$cov)[arrIndex.vcov[,"col"]]
        X.vcov <- rownames(object$cov)[arrIndex.vcov[,"row"]]
    
        name.vcov <- paste0(Y.vcov,
                            lava.options()$symbols[2],
                            X.vcov)
    
        value.vcov <- object$covfix[index.vcov]
        param.vcov <- unlist(Map(function(iPar,iFix,iName){if(iFix){NA}else if(!is.na(iPar)){iPar} else {iName}},
                                 iPar = object$covpar[index.vcov],
                                 iFix = !is.na(value.vcov),
                                 iName = name.vcov)
                             )

        index.variance <- which(arrIndex.vcov[,1]==arrIndex.vcov[,2])
        ls.name$variance <- name.vcov[index.variance]
        n.variance <- length(ls.name$variance)
        ls.Y$variance <- Y.vcov[index.variance]
        ls.X$variance <- X.vcov[index.variance]
        ls.type$variance <- rep("variance", n.variance)
        ls.value$variance <- value.vcov[index.variance]
        ls.param$variance <- param.vcov[index.variance]
        ls.marginal$variance <- ls.name$variance %in% paste0(exogenous(object),lava.options()$symbols[2],exogenous(object))
    
        index.covariance <- which(arrIndex.vcov[,1]!=arrIndex.vcov[,2])
        ls.name$covariance <- name.vcov[index.covariance]
        n.covariance <- length(ls.name$covariance)
        ls.Y$covariance <- Y.vcov[index.covariance]
        ls.X$covariance <- X.vcov[index.covariance]
        ls.type$covariance <- rep("covariance", n.covariance)
        ls.value$covariance <- value.vcov[index.covariance]
        ls.param$covariance <- param.vcov[index.covariance]
        ls.marginal$covariance <- rep(FALSE, n.covariance)
    }
    
    ## ** external coefficients
    n.external <- length(object$expar)
    if(n.external>0){
        ls.name$external <- names(object$expar)
        ls.type$external <- rep("external", n.external)

        ls.X$external <- rep(NA,n.external)
        for(iX in names(object$attributes$ordinalparname)){ ## iX <- "X1"
            ls.X$external[ls.name$external %in% object$attributes$ordinalparname[[iX]]] <- iX
        }
        ls.Y$external <- rep(NA,n.external)
        ls.value$external <- unlist(object$exfix)
        ls.param$external <- unlist(Map(function(iPar,iFix,iName){if(iFix){NA}else if(!is.na(iPar)){iPar} else {iName}},
                                        iPar = rep(NA,n.external),
                                        iFix = !is.na(ls.value$external),
                                        iName = ls.name$external)
                                    )
        ls.marginal$external <-  rep(FALSE, n.external)
    }

    ## ** merge
    df.param <- data.frame(name = unlist(ls.name),
                           Y = unlist(ls.Y),
                           X = unlist(ls.X),
                           data = unlist(ls.X),
                           type = unlist(ls.type),
                           value = unlist(ls.value),
                           param = unlist(ls.param),
                           marginal = unlist(ls.marginal),
                           stringsAsFactors = FALSE)
    df.param[df.param$X %in% latent(object),"data"] <- NA
    
    ## ** categorical variables
    if(!is.null(object$attributes$ordinalparname)){
        resCar <- defineCategoricalLink(object, link = df.param$name, data = data)
        
        ## normal coefficients
        resCar.Nexternal <- subset(resCar,
                                   subset = is.na(externalLink),
                                   select = c("link","type","factitious","level","originalLink","externalLink"))
        ## rename according to the main data frame
        match.tempo <- match(c("link","type"),
                             names(resCar.Nexternal))
        names(resCar.Nexternal)[match.tempo] <- c("name","distribution") 
        
        df.Nexternal <- merge(subset(df.param, subset = type != "external"),
                              resCar.Nexternal, by = "name")

        ## external coefficients
        resCar.external <- subset(resCar,
                                  subset = !is.na(externalLink),
                                  select = c("link", "endogenous", "exogenous", "type", "factitious", "level", "originalLink", "externalLink"))
        resCar.external$X <- paste0(resCar.external$exogenous,
                                    resCar.external$level)

        ## rename according to the main data frame
        match.tempo <- match(c("link","endogenous","exogenous","type"),
                             names(resCar.external))
        names(resCar.external)[match.tempo] <- c("name","Y","data","distribution")
        resCar.external$param <- resCar.external$name

        for(iCol in c("type","value","marginal")){ # iCol <- "type"
            name2col <- stats::setNames(df.param[[iCol]],df.param$name)
            resCar.external[,iCol] <- name2col[resCar.external$originalLink]
        }
        df.param <- rbind(resCar.external[,names(df.Nexternal),drop=FALSE],
                          df.Nexternal)
    }else{
        df.param$factitious <- FALSE
        df.param$level <- as.character(NA)
        df.param$externalLink <-  as.character(NA)
        df.param$originalLink <- df.param$name
    }

    ## ** merge with lava
    coef.lava <- coef(object)
    name.coef <- names(coef.lava)

    index.keep <- which(df.param$type!="external" & df.param$factitious == FALSE & df.param$marginal == FALSE)
    df.param$detail <- as.character(NA)
    df.param[index.keep, "detail"] <- detailName(object,
                                                 name.coef = df.param[index.keep, "name"],
                                                 type.coef = df.param[index.keep, "type"])
    df.param$lava <- name.coef[match(df.param$originalLink,coef.lava)]
    df.param <- df.param[order(df.param$type,df.param$detail,df.param$name),,drop=FALSE]
    rownames(df.param) <- NULL

    ## ** export
    if(as.lava){
        ## add extra mean as links
        vec.extra <- unique(stats::na.omit(df.param$externalLink))
        if(length(vec.extra)>0){
            df.extra <- data.frame(name = vec.extra, type = "extra",
                                   lava = name.coef[match(vec.extra,coef.lava)],
                                   stringsAsFactors = FALSE)
            df.param <- rbind(df.param[,c("name", "type", "lava")],
                              df.extra)
        }
        
        ## 
        out <- subset(df.param, subset = !is.na(lava), select = c("type", "name"))
        out <- stats::setNames(out$type, out$name)
        out <- out[!duplicated(names(out))]
        return(out[coef.lava])    
    }else{
        df.param$detail <- factor(df.param$detail,
                                  levels = c("nu","alpha","K","Gamma","Lambda","B","Sigma_var","Sigma_cov","sigma2","sigma2k","cor","Psi_var","Psi_cov",NA))
        df.param <- df.param[order(df.param$detail,df.param$param),]
        if(attr.type){
            iTempo <- detailType(df.param, endogenous = endogenous(object), latent = latent(object))
            attr(df.param, "byType") <- iTempo$byType
            attr(df.param, "byTypebyEndo") <- iTempo$byTypebyEndo
        }
        return(df.param)
    }
}

## * coefType.lvmfit
#' @rdname coefType
#' @export
coefType.lvmfit <- function(object, as.lava = TRUE, ...){ 

    ## ** find type of the coefficients in the original model
    df.param <- coefType(object$model0, as.lava = FALSE, ...)
    
    ## ** export
    if(as.lava){
        out <- subset(df.param, subset = !is.na(lava), select = c("type", "name"))
        out <- stats::setNames(out$type, out$name)
        coef.lava <- names(stats::coef(object))
        return(out[coef.lava])    
    }else{        
        return(df.param)
    }
}

## * coefType.multigroup
#' @rdname coefType
#' @export
coefType.multigroup <- function(object, as.lava = TRUE, attr.type = FALSE, ...){

    n.model <- length(object$lvm)
    df.param <- NULL
    
    for(iModel in 1:n.model){ # iModel <- 2

        df.param <- rbind(df.param,
                          cbind(coefType(object$lvm[[iModel]], as.lava = FALSE), model = iModel)
                          )
      
    }

    df.param$name <- paste0(df.param$model,"@",df.param$name)
    
    ## *** export
    if(as.lava){        
        out <- subset(df.param, subset = !is.na(lava), select = c("type", "name"))
        out <- stats::setNames(out$type, out$name)
        return(out)
    }else{
        if(attr.type){
            iTempo <- detailType(df.param, endogenous = endogenous(object), latent = latent(object))
            attr(df.param, "byType") <- iTempo$byType
            attr(df.param, "byTypebyEndo") <- iTempo$byTypebyEndo
        }
        return(df.param)
    }
}


## * detailName (needed for coefType)
detailName <- function(object, name.coef, type.coef){
    ls.links <- initVarLinks(name.coef)
    index.loading <- setdiff(which(ls.links$var2 %in% latent(object)),
                             which(type.coef %in% c("covariance","variance")))
    if(length(index.loading)>0){
        type.coef[index.loading] <- "loading"
    }

    index.measurement <- which(ls.links$var1 %in% endogenous(object))
    if(length(index.measurement)>0){
        type.coef[index.measurement] <- as.character(factor(type.coef[index.measurement],
                                                       levels = c("intercept","regression","loading","covariance","variance"),
                                                       labels = c("nu","K","Lambda","Sigma_cov","Sigma_var")))
    }
    index.structural <- setdiff(1:length(type.coef),index.measurement)
    if(length(index.structural)>0){
        type.coef[index.structural] <- as.character(factor(type.coef[index.structural],
                                                      levels = c("intercept","regression","loading","covariance","variance"),
                                                      labels = c("alpha","Gamma","B","Psi_cov","Psi_var")))
    }
    return(type.coef)
}

## * detailType (needed for coefType)
detailType <- function(data, endogenous, latent){
    n.endogenous <- length(endogenous)
    n.latent <- length(latent)
    
    data <- data[!is.na(data$param),]
    Utype <- levels(droplevels(data$detail))
    
    out <- list()
    out$byType <- tapply(data$name,droplevels(data$detail),function(iVec){list(unique(iVec))})
    out$byTypebyEndo <- list()

    ## ** nu
    if("nu" %in% data$detail){
        data.nu <- data[data$detail=="nu",,drop=FALSE]
        data.nu <- data.nu[match(endogenous,data.nu$Y),]
        out$byTypebyEndo$nu <- setNames(data.nu$param,data.nu$Y)
    }

    ## ** Lambda
    if("Lambda" %in% data$detail){
        data.Lambda <- data[data$detail=="Lambda",,drop=FALSE]
        data.Lambda <- data.Lambda[match(endogenous,data.Lambda$Y),]
        out$byTypebyEndo$Lambda <- setNames(data.Lambda$param,data.Lambda$Y)
    }
    
    ## ** K and X
    if("K" %in% data$detail){
        data.K <- data[data$detail=="K",,drop=FALSE]
        data.K <- data.K[match(endogenous,data.K$Y),]
        out$byTypebyEndo$K <- setNames(lapply(endogenous, function(iE){data.K[data.K$Y == iE,"param"]}), endogenous)
        out$byTypebyEndo$XK <- setNames(lapply(endogenous, function(iE){data.K[data.K$Y == iE,"X"]}), endogenous)
    }

    ## ** Sigma
    if("Sigma_var" %in% data$detail || "Sigma_cov" %in% data$detail){
        data.Sigma <- data[data$detail %in% c("Sigma_var","Sigma_cov"),,drop=FALSE]
        out$byTypebyEndo$Sigma <- matrix(NA, nrow = n.endogenous, ncol = n.endogenous,
                                          dimnames = list(endogenous, endogenous))

        for(iSigma in 1:NROW(data.Sigma)){
            out$byTypebyEndo$Sigma[data.Sigma$X[iSigma],data.Sigma$Y[iSigma]] <- data.Sigma$param[iSigma]
            out$byTypebyEndo$Sigma[data.Sigma$Y[iSigma],data.Sigma$X[iSigma]] <- data.Sigma$param[iSigma]
        }
    }
    if("cor" %in% data$detail){
        data.cor <- data[data$detail %in% "cor",,drop=FALSE]
        out$byTypebyEndo$cor <- matrix(NA, nrow = n.endogenous, ncol = n.endogenous,
                                       dimnames = list(endogenous, endogenous))
        for(iCor in 1:NROW(data.cor)){
            out$byTypebyEndo$cor[data.cor$X[iCor],data.cor$Y[iCor]] <- data.cor$param[iCor]
            out$byTypebyEndo$cor[data.cor$Y[iCor],data.cor$X[iCor]] <- data.cor$param[iCor]
        }
    }
    if("sigma2" %in% data$detail){
        out$byTypebyEndo$sigma2 <- matrix("sigma2", nrow = n.endogenous, ncol = n.endogenous,
                                          dimnames = list(endogenous, endogenous))
        if(is.null(out$byTypebyEndo$cor)){
            out$byTypebyEndo$sigma2[upper.tri(out$byTypebyEndo$sigma2)] <- NA
            out$byTypebyEndo$sigma2[lower.tri(out$byTypebyEndo$sigma2)] <- NA
        }
    }
    if("sigma2k" %in% data$detail){
        data.sigma2k <- data[data$detail %in% "sigma2k",,drop=FALSE]

        out$byTypebyEndo$sigma2k <- setNames(lapply(1:NROW(data.sigma2k), function(iP){
            matrix(NA, nrow = n.endogenous, ncol = n.endogenous,
                   dimnames = list(endogenous, endogenous))
        }), data.sigma2k$param)

        
        for(iSigma in 1:NROW(data.sigma2k)){ # iSigma <- 1
            out$byTypebyEndo$sigma2k[[iSigma]][data.sigma2k$Y[iSigma],] <- data.sigma2k$param[iSigma]
            out$byTypebyEndo$sigma2k[[iSigma]][,data.sigma2k$Y[iSigma]] <- data.sigma2k$param[iSigma]
        }
    }

    ## ** alpha
    if("alpha" %in% data$detail){
        data.alpha <- data[data$detail=="alpha",,drop=FALSE]
        data.alpha <- data.alpha[match(latent,data.alpha$Y),]
        out$byTypebyEndo$alpha <- setNames(data.alpha$param,data.alpha$Y)
    }

    ## ** B
    if("B" %in% data$detail){
        data.B <- data[data$detail %in% c("B"),,drop=FALSE]
        out$byTypebyEndo$B <- matrix(NA, nrow = n.latent, ncol = n.latent,
                                     dimnames = list(latent, latent))
        for(iB in 1:NROW(data.B)){ # iSigma <- 1
            out$byTypebyEndo$B[data.B$X[iB],data.B$Y[iB]] <- data.B$param[iB]
        }
    
    }
    
    ## ** Gamma and X
    if("Gamma" %in% data$detail){
        data.Gamma <- data[data$detail=="Gamma",,drop=FALSE]
        data.Gamma <- data.Gamma[match(latent,data.Gamma$Y),]
        out$byTypebyEndo$Gamma <- setNames(lapply(latent, function(iE){data.Gamma[data.Gamma$Y == iE,"param"]}), latent)
        out$byTypebyEndo$XGamma <- setNames(lapply(latent, function(iE){data.Gamma[data.Gamma$Y == iE,"X"]}), latent)
    }

    ## ** Psi
    if("Psi_var" %in% data$detail || "Psi_cov" %in% data$detail){
        data.Psi <- data[data$detail %in% c("Psi_var","Psi_cov"),,drop=FALSE]
        out$byTypebyEndo$Psi <- matrix(NA, nrow = n.latent, ncol = n.latent,
                                         dimnames = list(latent, latent))

        for(iPsi in 1:NROW(data.Psi)){
            out$byTypebyEndo$Psi[data.Psi$X[iPsi],data.Psi$Y[iPsi]] <- data.Psi$param[iPsi]
            out$byTypebyEndo$Psi[data.Psi$Y[iPsi],data.Psi$X[iPsi]] <- data.Psi$param[iPsi]
        }
    }
        
    return(out)
}



##----------------------------------------------------------------------
### coefType.R ends here

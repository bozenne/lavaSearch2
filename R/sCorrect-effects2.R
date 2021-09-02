### effects2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  4 2019 (10:28) 
## Version: 
## Last-Updated: sep  2 2021 (14:14) 
##           By: Brice Ozenne
##     Update #: 229
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * effects2 (documentation)
#' @title Effects From a Fitted Model After Small Sample Correction 
#' @description Test whether a path in the latent variable model correspond to a null effect.
#' Similar to \code{lava::effects} but with small sample correction (if any).
#' So far it only work for paths composed of two edges.
#' @name effects2
#'
#' @param object an object that inherits from lvmfit.
#' @param link [character vector] The path for which the effect should be assessed (e.g. \code{"A~B"}),
#' i.e. the effect of the right variable (B) on the left variable (A). 
#' @param conf.level [numeric, 0-1] confidence level of the interval.
#' @param df [logical] should the degree of freedoms of the Wald statistic be computed using the Satterthwaite correction?
#' Otherwise the degree of freedoms are set to \code{Inf}, i.e. a normal distribution is used instead of a Student's t distribution when computing the p-values.
#' @param ssc [logical] should the standard errors of the coefficients be corrected for small sample bias? Argument passed to \code{sCorrect}.
#' @param ... [logical] arguments passed to lower level methods.
#' 
#' @concept small sample inference
#' @export
`effects2` <-
  function(object, ...) UseMethod("effects2")

## * effects2 (examples)
## TODO

## * effects2.lvmfit
#' @rdname effects2
#' @export
effects2.lvmfit <- function(object, ssc = lava.options()$ssc, df = lava.options()$df, ...){

    object.SSC <- sCorrect(object, ssc = ssc, df = df)    
    return(effects2(object.SSC, ...))

}

## * effects2.sCorrect
#' @rdname effects2
#' @export
effects2.sCorrect <- function(object, link, conf.level = 0.95, robust = FALSE, ...){

    ## ** extract information
    allCoef.type <- coefType(object, as.lava = FALSE)
    n.hypo <- length(link)
    name.hypo <- names(link)
    if(is.null(name.hypo)){name.hypo <- link}
    test.df <- identical(object$sCorrect$df,"Satterthwaite")
    
    object.coef <- coef2(object)
    object.iid <- iid2(object, robust = robust)
    object.vcov.param <- vcov2(object)
    if(test.df){
        object.dVcov.param <- object$sCorrect$dVcov.param
    }
        
    name.coef <- names(object.coef)
    n.coef <- length(name.coef)
    n.obs <- NROW(object.iid)

    ## ** extract coefficient and contrast matrix
    null <- rep(NA, n.hypo)
    ls.contrast <- vector(mode = "list", length = n.hypo)
    for(iH in 1:n.hypo){ # iH <- 1
        iTempo.eq <- strsplit(link[iH], split = "=", fixed = TRUE)[[1]]
        if(length(iTempo.eq)==1){ ## set null to 0 when second side of the equation is missing
            iTempo.eq <- c(iTempo.eq,"0")
        }

        null[iH] <- as.numeric(trim(iTempo.eq[2]))
        iRh.plus <- strsplit(iTempo.eq[[1]], split = "+", fixed = TRUE)[[1]]
        iRh <- trim(unlist(sapply(iRh.plus, strsplit, split = "-", fixed = TRUE)))
        iRh <- iRh[iRh!="",drop=FALSE]
                            
        ls.iRh <- lapply(strsplit(iRh, split = "*", fixed = TRUE), trim)
        iN.tempo <- length(ls.iRh)
        ls.contrast[[iH]] <- setNames(rep(0, iN.tempo), sapply(ls.iRh, function(x){tail(x,1)}))
        for(iCoef in 1:iN.tempo){ # iCoef <- 2
            if(length(ls.iRh[[iCoef]])==1){
                iFactor <- 1
                iName <- ls.iRh[[iCoef]][1]                
            }else{
                iFactor <- as.numeric(ls.iRh[[iCoef]][1])
                iName <- ls.iRh[[iCoef]][2]
            }

            ## identify if it is a minus sign
            iBeforeCoef <- strsplit(iTempo.eq[[1]], split = ls.iRh[iCoef])[[1]][1]
            if(iCoef > 1){
                iBeforeCoef <- strsplit(iBeforeCoef, split = ls.iRh[iCoef-1])[[1]][2]
            }
            test.sign <- length(grep("-",iBeforeCoef))>0
            ls.contrast[[iH]][iName] <- c(1,-1)[test.sign+1] * iFactor
        }
                    
    }

    name.link <- unique(unlist(lapply(ls.contrast,names)))
    n.link <- length(name.link)
    contrast <- matrix(0, nrow = n.hypo, ncol = n.link,
                       dimnames = list(name.hypo, name.link))
    for(iH in 1:n.hypo){ # iH <- 1
        contrast[iH, names(ls.contrast[[iH]])] <- unname(ls.contrast[[iH]])
    }

    ## ** identify paths
    ls.link <- setNames(vector(mode = "list", length = n.link), name.link)
    
    for(iL in 1:n.link){ ## iL <- 1
        iLink <- name.link[iL]
        iPath <- lava::path(object, to = as.formula(iLink))

        if(length(iPath$path[[1]])==0){
            stop("No path found \n")
        }else{
            iNode <- iPath$path[[1]]
            iN.node <- length(iNode)
            ls.link[[iL]] <- paste0(iNode[-1], lava.options()$symbols[1], iNode[-iN.node], collpase = "")

            if(any(ls.link[[iL]] %in% allCoef.type$name == FALSE)){
                stop("Part of the path could not be identified \n")
            }
            if(any(allCoef.type$type[allCoef.type$name %in% ls.link[[iL]]] != "regression")){
                stop("Part of the path does not correspond to a regression link \n")
            }
        }
    }
    ## ** apply chain rule over the path
    vec.beta <- setNames(rep(NA, length = n.link), name.link)
    vec.sd <- setNames(rep(NA, length = n.link), name.link)
    vec.df <- setNames(rep(NA, length = n.link), name.link)
    M.iid <- matrix(NA, nrow = n.obs, ncol = n.link,
                    dimnames = list(NULL, name.link))
    if(test.df){
        M.dVcov <- matrix(NA, nrow = n.coef, ncol = n.link,
                          dimnames = list(name.coef, name.link))
    }
    
    for(iL in 1:n.link){ ## iL <- 1
        iLink <- ls.link[[iL]]
        vec.beta[iL] <- object.coef[iLink[1]]
        M.iid[,iL] <- object.iid[,iLink[1]]
        if(test.df){
            M.dVcov[,iL] <- object.dVcov.param[iLink[1],iLink[1],]
        }
        
        if(length(iLink)>1){
            for(iL2 in 2:length(iLink)){ ## iL2 <- 2
                iType <- allCoef.type[which(iLink[iL2] == allCoef.type$name),]
                if(!is.na(iType$value)){
                    vec.beta[iL] <- vec.beta[iL] * iType$value
                    M.iid[,iL] <- M.iid[,iL] * iType$value
                    if(test.df){
                        M.dVcov[,iL] <- M.dVcov[,iL] * iType$value^2
                    }
                }else{
                    betaOLD <- vec.beta[iL]
                    betaNEW <- object.coef[iLink[iL2]]
                    vec.beta[iL] <- betaOLD * betaNEW

                    iidOLD <- M.iid[,iL]
                    iidNEW <- object.iid[,iLink[iL2]]
                    M.iid[,iL] <- iidOLD * betaNEW + iidNEW * betaOLD

                    if(test.df){
                        varOLD <- sum(iidOLD^2)
                        varNEW <- sum(iidNEW^2)
                        covOLDNEW <- sum(iidOLD*iidNEW)

                        dVcovOLD <- M.dVcov[,iL]
                        dVcovNEW <- object.dVcov.param[iLink[iL2],iLink[iL2],]
                        dVcovCOV <- object.dVcov.param[iLink[iL],iLink[iL2],]
                        object.dVcov.param[iLink[iL],iLink[iL2],]
                        object.dVcov.param[iLink[iL2],iLink[iL],]
                    
                        browser()

                        Ilink1 <- as.numeric(keep.param %in% link1)
                        Ilink2 <- as.numeric(keep.param %in% link2)

                        ## M.var <- varOLD * betaNEW^2 + varNEW * betaOLD^2 + 2 * covOLDNEW * betaNEW * betaOLD
                        ## so dM.var/dtheta =
                        M.dVcov[,iL] <- (dVcovOLD * betaNEW^2 + 2 * varOLD * betaNEW) + (dVcovNEW * betaOLD^2 + 2 * varNEW * betaOLD) + 2 * dVcovOLDNEW * 
                            dSigma1 <- M.dVcov[,iL]
                        dSigma2 <- object.dVcov.param[iLink[iL2],iLink[iL2],]


                        ## effect.var <- as.double(Sigma[link1,link1] * mu[link2]^2 + Sigma[link2,link2] * mu[link1]^2 + 2 * Sigma[link1,link2] * mu[link1] * mu[link2])
                        ## dvar1 <- dSigma[link1,link1,] * mu[link2]^2 + Sigma[link1,link1] * 2 * Ilink2 * mu[link2]
                        ## dvar2 <- dSigma[link2,link2,] * mu[link1]^2 + Sigma[link2,link2] * 2 * Ilink1 * mu[link1]
                        ## dvar12 <- 2 * dSigma[link1,link2,] * mu[link1] * mu[link2] + 2 * Sigma[link1,link2] * (Ilink2 * mu[link2] + mu[link1] * Ilink1)
                        ## dvar <- dvar1 + dvar2 + dvar12
                    }
                    

                }
            }
        }

    }

    ## ** compute df
    if(test.df){
        vec.df <- 2*vec.var^2 / rowSums(t(M.dVcov) %*% object.vcov.param * t(M.dVcov))
    }

    ## ** gather everything in glht object
    out <- list(model = object,
                linfct = contrast,
                rhs = null,
                coef = vec.beta,
                vcov = crossprod(M.iid),
                df = if(test.df){vec.df}else{0},
                alternative = "two.sided",
                type = NULL,
                robust = robust,
                ssc = object$sCorrect$ssc$type,
                global = NULL)
    class(out) <- c("glht2","glht")

    ## ** export
    return(out)

}

## ## * .deltaMethod_product
## .deltaMethod_product <- function(mu,Sigma,dSigma,link){
##     link1 <- link[1]
##     link2 <- link[2]

##     effect <- as.double(prod(mu[link]))
##     effect.var <- as.double(Sigma[link1,link1] * mu[link2]^2 + Sigma[link2,link2] * mu[link1]^2 + 2 * Sigma[link1,link2] * mu[link1] * mu[link2])
##     effect.se <- sqrt(effect.var)
##     effect.Wald <- effect/effect.se

##     if(!is.null(dSigma)){
##         keep.param <- dimnames(dSigma)[[3]]
##         Ilink1 <- as.numeric(keep.param %in% link1)
##         Ilink2 <- as.numeric(keep.param %in% link2)
    
##         dvar1 <- dSigma[link1,link1,] * mu[link2]^2 + Sigma[link1,link1] * 2 * Ilink2 * mu[link2]
##         dvar2 <- dSigma[link2,link2,] * mu[link1]^2 + Sigma[link2,link2] * 2 * Ilink1 * mu[link1]
##         dvar12 <- 2 * dSigma[link1,link2,] * mu[link1] * mu[link2] + 2 * Sigma[link1,link2] * (Ilink2 * mu[link2] + mu[link1] * Ilink1)
##         dvar <- dvar1 + dvar2 + dvar12

##         effect.df <- 2 * effect.var^2 / (t(dvar) %*% Sigma[keep.param,keep.param,drop=FALSE] %*% dvar)[1,1]
##     }else{
##         effect.df <- Inf
##     }

##     ## ** export
##     iOut <- c("estimate" = as.double(effect),
##               "std.error" = as.double(effect.se),
##               "df" = as.double(effect.df),
##               "ci.lower" = as.numeric(NA),
##               "ci.upper" = as.numeric(NA),
##               "statistic" = as.double(effect.Wald),
##               "p.value" = as.numeric(NA)
##               )
##     if(is.infinite(iOut["df"])){                       
##         iOut["p.value"] <- as.double(2*(1-pnorm(abs(iOut["statistic"]))))
##     }else{
##         iOut["p.value"] <- as.double(2*(1-pt(abs(iOut["statistic"]), df = iOut["df"])))
##     }
##     return(iOut)
## }

######################################################################
### effects2.R ends here

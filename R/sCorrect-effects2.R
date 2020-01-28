### effects2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  4 2019 (10:28) 
## Version: 
## Last-Updated: jan 27 2020 (16:16) 
##           By: Brice Ozenne
##     Update #: 143
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
#' @param transform [function] function to backtransform the estimates and the associated confidence intervals.
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
    return(effects(object.SSC, ...))

}

## * effects2.sCorrect
#' @rdname effects2
effects2.sCorrect <- function(object, link, conf.level = 0.95, transform = NULL, ...){

    ## ** compute product
    n.link <- length(link)

    name.coef <- names(coef(object))
    link.direct <- link[link %in% name.coef]
    link.other <- setdiff(link, link.direct)

    object.summary2 <- summary2(object, ssc = ssc, df = df)$coef 

    out <- NULL
    if(length(link.direct)>0){        
        out <- rbind(out,
                     object.summary2[link.direct,])
    }

    if(length(link.other)>0){
        allCoef <- coefType(object, as.lava = FALSE)
        mu <- setNames(object.summary2[,"Estimate"],rownames(object.summary2))
        mu.se <- setNames(object.summary2[,"Std. Error"],rownames(object.summary2))
        mu.df <- setNames(object.summary2[,"df"],rownames(object.summary2))
        Sigma <- vcov2(object)
        dSigma <- object$sCorrect$dVcov.param
        
        for(iL in 1:length(link.other)){ ## iL <- 1
            iLink <- link.other[iL]
            iPath <- lava::path(object, to = as.formula(iLink))
            if(length(iPath$path[[1]])==0){
                stop("No path found \n")
            }else{
                iNode <- iPath$path[[1]]
                iN.node <- length(iNode)
                iLink <- paste0(iNode[-1], lava.options()$symbols[1], iNode[-iN.node], collpase = "")

                if(any(iLink %in% allCoef$name == FALSE)){
                    stop("Part of the path could not be identified \n")
                }
                if(any(allCoef$type[allCoef$name %in% iLink] != "regression")){
                    stop("Part of the path does not correspond to a regression link \n")
                }
                if(any(length(iLink)!=2)){
                    stop("Only implemented for path of length 2 \n")
                }                

                test.noconstrain <- iLink %in% allCoef$originalLink
                if(all(test.noconstrain)){
                    out <- rbind(out,
                                 .deltaMethod_product(mu = mu, Sigma = Sigma, dSigma = dSigma, link = iLink)
                                 )
                }else{
                    iLink.NA <- iLink[test.noconstrain==FALSE]
                    iLink.NNA <- iLink[test.noconstrain==TRUE]
                    iEffect <- prod(mu[iLink])
                    iEffect.se <- mu.se[iLink.NNA] * mu[iLink.NA]
                    iEffect.df <- mu.df[iLink.NNA]

                    iRow <- c("estimate" = as.double(iEffect),                                  
                              "std.error" = as.double(iEffect.se),
                              "df" = as.double(iEffect.df),
                              "ci.lower" = as.numeric(NA),
                              "ci.upper" = as.numeric(NA),
                              "statistic" = as.double(iEffect/iEffect.se),
                              "p.value" = as.numeric(NA)
                              )

                    if(is.infinite(iRow["df"])){                       
                        iRow["p.value"] <- as.double(2*(1-pnorm(abs(iRow["statistic"]))))
                    }else{
                        iRow["p.value"] <- as.double(2*(1-pnorm(abs(iRow["statistic"]), df = iRow["df"])))
                    }
                    out <- rbind(out,iRow)                    
                }
                
            }

        }            
        
    }
    rownames(out) <- c(link.direct,link.other)

    ## ** transform
    if(is.null(transform)){
        transform <- function(x){x}
        dtransform <- function(x){x}
    }else if(!is.null(attr(transform,"derivative"))){
        dtransform <- attr(transform,"derivative")
    }else{
        dtransform <- function(x){numDeriv::jacobian(transform, x)[1,1]}
    }

    if(all(is.infinite(out[,"df"]))){
        out[,"ci.lower"] <- sapply(out[,"estimate"] + qnorm((1-conf.level)/2) * out[,"std.error"], transform)
        out[,"ci.upper"] <- sapply(out[,"estimate"] + qnorm(conf.level + (1-conf.level)/2) * out[,"std.error"], transform)
    }else{
        out[,"ci.lower"] <- sapply(out[,"estimate"] + qt((1-conf.level)/2, df = out[,"df"]) * out[,"std.error"], transform)
        out[,"ci.upper"] <- sapply(out[,"estimate"] + qt(conf.level + (1-conf.level)/2, df = out[,"df"]) * out[,"std.error"], transform)
    }
    out[,"std.error"] <- out[,"std.error"]*dtransform(out[,"estimate"])
    out[,"estimate"] <- transform(out[,"estimate"])
    
    ## ** export    
    return(out)

}

## * effects.sCorrect
#' @rdname effects2
effects.sCorrect <- effects2.sCorrect

## * .deltaMethod_product
.deltaMethod_product <- function(mu,Sigma,dSigma,link){
    link1 <- link[1]
    link2 <- link[2]

    effect <- as.double(prod(mu[link]))
    effect.var <- as.double(Sigma[link1,link1] * mu[link2]^2 + Sigma[link2,link2] * mu[link1]^2 + 2 * Sigma[link1,link2] * mu[link1] * mu[link2])
    effect.se <- sqrt(effect.var)
    effect.Wald <- effect/effect.se

    if(!is.null(dSigma)){
        keep.param <- dimnames(dSigma)[[3]]
        Ilink1 <- as.numeric(keep.param %in% link1)
        Ilink2 <- as.numeric(keep.param %in% link2)
    
        dvar1 <- dSigma[link1,link1,] * mu[link2]^2 + Sigma[link1,link1] * 2 * Ilink2 * mu[link2]
        dvar2 <- dSigma[link2,link2,] * mu[link1]^2 + Sigma[link2,link2] * 2 * Ilink1 * mu[link1]
        dvar12 <- 2 * dSigma[link1,link2,] * mu[link1] * mu[link2] + 2 * Sigma[link1,link2] * (Ilink2 * mu[link2] + mu[link1] * Ilink1)
        dvar <- dvar1 + dvar2 + dvar12

        effect.df <- 2 * effect.var^2 / (t(dvar) %*% Sigma[keep.param,keep.param,drop=FALSE] %*% dvar)[1,1]
    }else{
        effect.df <- Inf
    }

    ## ** export
    iOut <- c("estimate" = as.double(effect),
              "std.error" = as.double(effect.se),
              "df" = as.double(effect.df),
              "ci.lower" = as.numeric(NA),
              "ci.upper" = as.numeric(NA),
              "statistic" = as.double(effect.Wald),
              "p.value" = as.numeric(NA)
              )
    if(is.infinite(iOut["df"])){                       
        iOut["p.value"] <- as.double(2*(1-pnorm(abs(iOut["statistic"]))))
    }else{
        iOut["p.value"] <- as.double(2*(1-pnorm(abs(iOut["statistic"]), df = iOut["df"])))
    }
    return(iOut)
}

######################################################################
### effects2.R ends here

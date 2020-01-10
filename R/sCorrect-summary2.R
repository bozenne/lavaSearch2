### sCorrect-summary2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 10 2017 (10:57) 
## Version: 
## Last-Updated: jan  8 2020 (16:38) 
##           By: Brice Ozenne
##     Update #: 348
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - summary2
#' @title Summary After Small Sample Correction
#' @description Summary of the estimates, standard errors, and p-values based on Wald tests for \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} objects.
#' Similar to \code{summary} but with small sample correction (if any).
#' @name summary2
#'
#' @param object a \code{lm}, \code{gls}, \code{lme} or \code{lvm} object.
#' @param digit [integer > 0] the number of decimal places to use when displaying the summary.
#' @param df [logical] should the degree of freedoms of the Wald statistic be computed using the Satterthwaite correction?
#' Otherwise the degree of freedoms are set to \code{Inf}, i.e. a normal distribution is used instead of a Student's t distribution when computing the p-values.
#' @param bias.correct [logical] should the standard errors of the coefficients be corrected for small sample bias?
#' See \code{\link{sCorrect}} for more details.
#' @param robust [logical] should the robust standard errors be used instead of the model based standard errors?
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param ... arguments passed to the \code{summary} method of the object.
#' 
#' @seealso \code{\link{sCorrect}} for more detail about the small sample correction.
#'
#' @details \code{summary2} is the same as \code{summary}
#' except that it first computes the small sample correction (but does not store it).
#' So if \code{summary2} is to be called several times,
#' it is more efficient to pre-compute the quantities for the small sample correction
#' using \code{sCorrect} and then call \code{summary2}.
#' 
#' @examples
#' m <- lvm(Y~X1+X2)
#' set.seed(10)
#' d <- lava::sim(m, 2e1)
#'
#' ## Gold standard
#' summary(lm(Y~X1+X2, d))$coef
#' 
#' ## gls models
#' library(nlme)
#' e.gls <- gls(Y~X1+X2, data = d, method = "ML")
#' summary(e.gls)$tTable
#' sCorrect(e.gls, cluster = 1:NROW(d)) <- FALSE ## no small sample correction
#' summary2(e.gls)$tTable
#' 
#' sCorrect(e.gls, cluster = 1:NROW(d)) <- TRUE ## small sample correction
#' summary2(e.gls)$tTable
#' 
#' ## lvm models
#' e.lvm <- estimate(m, data = d)
#' summary(e.lvm)$coef
#' 
#' sCorrect(e.lvm) <- FALSE ## no small sample correction
#' summary2(e.lvm)$coef
#' 
#' sCorrect(e.lvm) <- TRUE ## small sample correction
#' summary2(e.lvm)$coef
#' 
#' @concept small sample inference
#' @export
`summary2` <-
  function(object, ssc, df, digit, robust, ...) UseMethod("summary2")

## * summary2.lm
#' @rdname summary2
#' @method summary2 lm
#' @export
summary2.lm <- function(object, ssc = lava.options()$ssc, df = lava.options()$df,
                        digit = max(3, getOption("digit")),
                        robust = FALSE, ...){
    if(is.null(object$sCorrect) || !identical(object$sCorrect$ssc$type, ssc) || !identical(object$sCorrect$df, df)){
        object <- sCorrect(object, ssc = ssc, df = df)
    }
    ## ** perform Wald test
    name.param <- names(coef2(object, ssc = NA))
    n.param <- length(name.param)

    tTable.all <- compare2(object,
                           linfct = name.param,
                           robust = robust,
                           ssc = ssc, df = df,
                           F.test = FALSE,
                           as.lava = FALSE)
    tTable <- tTable.all[1:n.param,c("estimate","std","statistic","p-value","df")]
    dimnames(tTable) <- list(name.param,
                             c("Value","Std.Error","t-value","p-value","df")
                             )

    ## ** get summary
    class(object) <- setdiff(class(object),c("sCorrect"))
    object.summary <- summary(object, digits = digit, ...)
    
    ## ** update summary
    object.summary$coefficients <- tTable
    object.summary$residuals <- object$sCorrect$residuals[,1]
    object.summary$sigma <- sqrt(tTable["sigma2","Value"])

    object.summary$cov.unscaled <- NULL
    object.summary$fstatistic <- NULL

    ## ** export
    return(object.summary)

}
## * summary2.gls
#' @rdname summary2
#' @method summary2 gls
#' @export
summary2.gls <- function(object, ssc = lava.options()$ssc, df = lava.options()$df,
                         digit = max(3, getOption("digit")),
                         robust = FALSE, ...){
    if(is.null(object$sCorrect) || !identical(object$sCorrect$ssc$type, ssc) || !identical(object$sCorrect$df, df)){
        object <- sCorrect(object, ssc = ssc, df = df)
    }

    ## ** perform Wald test
    name.param <- names(coef2(object, ssc = ssc))
    n.param <- length(name.param)

    tTable.all <- compare2(object,
                           linfct = name.param,
                           robust = robust,
                           ssc = ssc, df = df,
                           F.test = FALSE,
                           as.lava = FALSE)
    tTable <- tTable.all[1:n.param,c("estimate","std","statistic","p-value","df")]
    dimnames(tTable) <- list(name.param,
                             c("Value","Std.Error","t-value","p-value","df")
                             )

    ### ** get summary
    class(object) <- setdiff(class(object),c("sCorrect"))
    object.summary <- summary(object, digits = digit, ...)
    
    ### ** update summary
    object.summary$tTable <- tTable
    browser()
    object.summary$residuals <- quantile(residuals2(object, ssc = ssc, type = "normalized"), na.rm = TRUE)
    quantile(residuals(object), na.rm = TRUE)
    
    object.summary$sigma <- tTable["sigma2","Value"]
    browser()
    
    ### ** export
    return(object.summary)
}
## * summary2.lme
#' @rdname summary2
#' @method summary2 lme
#' @export
summary2.lme <- summary2.gls

## * summary2.lvmfit
#' @rdname summary2
#' @method summary2 lvmfit
#' @export
summary2.lvmfit <- function(object, ssc = lava.options()$ssc, df = lava.options()$df,
                            digit = max(3, getOption("digit")),
                            robust = FALSE, ...){
    if(is.null(object$sCorrect) || !identical(object$sCorrect$ssc$type, ssc) || !identical(object$sCorrect$df, df)){
        object <- sCorrect(object, ssc = ssc, df = df)
    }

    ## ** perform Wald test
    param <- coef2(object, ssc = ssc)
    name.param <- names(param)
    n.param <- length(param)

    table.all <- compare2(object,
                          linfct = name.param,
                          robust = robust,
                          cluster = cluster,
                          ssc = ssc, df = df,
                          F.test = FALSE,
                          as.lava = FALSE)
    table.coef <- table.all[1:n.param,c("estimate","std","statistic","p-value","df")]
    dimnames(table.coef) <- list(name.param,
                                 c("Estimate", "Std. Error", "t-value", "P-value", "df")
                                 )

    ## ** get summary
    class(object) <- setdiff(class(object),"sCorrect")
    object.summary <- summary(object, ...)
    if(!is.null(object$cluster) || inherits(object,"lvm.missing")){
        
        ## if(robust == FALSE){
        ##     stop("Can only print summary for robust standard errors \n",
        ##          "when the object contain a cluster variable \n")
        ## }
        colnames(object.summary$coef) <- c("Estimate","Std. Error","Z-value","P-value")
        object.summary$coef[,"Z-value"] <- NA

        colnames(object.summary$coefmat) <- c("Estimate","Std. Error","Z-value","P-value", "std.xy")
        object.summary$coefmat[,"Z-value"] <- ""
        
    }
    
    ## find digit
    vec.char <- setdiff(object.summary$coefmat[,"Estimate"],"")
    digit <- max(c(nchar(gsub(".","",vec.char,fixed = TRUE)))-1,1)

    ## ** update summary
    ## *** vcov
    object.summary$vcov <- attr(object$dVcov, "vcov.param")[name.param,name.param]    

    ## *** coef
    lava.rownames <- rownames(object.summary$coef)
    ## add rows corresponding to reference parameters
    missing.rows <- setdiff(lava.rownames,rownames(table.coef))
    if(length(missing.rows)>0){
        addon <- object.summary$coef[missing.rows,
                                     c("Estimate","Std. Error","Z-value","P-value"),
                                     drop=FALSE]
        colnames(addon)[3] <- "t-value"
        table.coef <- rbind(table.coef, cbind(addon,df=NA))
    }

    ## re-order table according to lava
    table.coef <- table.coef[intersect(lava.rownames,rownames(table.coef)),,drop=FALSE]
    ## remove unappropriate p.values
    lava.NApvalue <- which(is.na(object.summary$coef[,"P-value"]))
    table.coef[intersect(lava.rownames[lava.NApvalue],rownames(table.coef)),"P-value"] <- NA
    object.summary$coef <- table.coef
    
    
    ## *** coefmat
    name.label0 <- trimws(rownames(CoefMat(object, labels = 0, level = 9)), which = "both")
    index.titleVariance <- which(name.label0=="Residual Variances:")
    if(length(index.titleVariance)>0){
        ## rename variance parameters from Y to Y~~Y
        index.vcov <- (index.titleVariance+1):length(name.label0)
        index.var <- setdiff(index.vcov,grep("~~",name.label0,fixed=TRUE)) ## exclude covariance parameters that are already correctly named
        name.label0[index.var] <- paste0(name.label0[index.var],lava.options()$symbols[2],name.label0[index.var])
    }

    table.coefmat <- object.summary$coefmat
    colnames(table.coefmat)[3:5] <- c("t-value","P-value","df")
    
    ## mimic lava:::CoefMat (called by lava:::summary.lvmfit)    
    e2add <- format(round(table.coef[,"Estimate"], max(1, digit - 1)), digits = digit - 1)
    e2add <- gsub(" NA","",e2add)
    sd2add <- format(round(table.coef[,"Std. Error"], max(1, digit - 1)), digits = digit - 1)
    sd2add <- gsub(" NA","",sd2add)
    df2add <- as.character(round(table.coef[,"df"],2))    
    df2add[is.na(df2add)] <- ""
    t2add <- format(round(table.coef[,"t-value"], max(1, digit - 1)), digits = digit - 1)
    t2add <- gsub(" NA","",t2add)

    p2add <- formatC(table.coef[,"P-value"], digits = digit - 1, format = "g",  preserve.width = "common", flag = "")
    p2add <- gsub(" NA","",p2add)
    p2add[table.coef[,"P-value"] < 1e-12] <- "  <1e-12"

    M2add <- cbind(e2add,sd2add,t2add,p2add,df2add)
    table.coefmat[match(rownames(table.coef), name.label0),] <- M2add

    table.coefmat[object.summary$coefmat[,"P-value"]=="","P-value"] <- ""
    object.summary$coefmat <- table.coefmat

    ## ** Export
    if(robust){
        colnames(object.summary$coefmat)[2] <- "robust SE"
        colnames(object.summary$coef)[2] <- "robust SE"
    }
    return(object.summary)    

}

## * summary2.sCorrect
#' @rdname summary2
summary2.sCorrect <- function(object, ssc = object$sCorrect$ssc$type, df = object$sCorrect$df,
                              digit = max(3, getOption("digit")),
                              robust = FALSE, ...){
    class(object) <- setdiff(class(object),"sCorrect")
    return(summary2(object, ssc = ssc, df = df, digit = digit, robust = robust, ...))
}

## * summary.sCorrect
#' @rdname summary2
summary.sCorrect <- summary2.sCorrect



##----------------------------------------------------------------------
### Scorrect-summary2.R ends here


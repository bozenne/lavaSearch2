
### sCorrect-summary2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 10 2017 (10:57) 
## Version: 
## Last-Updated: feb 11 2020 (17:24) 
##           By: Brice Ozenne
##     Update #: 456
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
#' @param robust [logical] should the robust standard errors be used instead of the model based standard errors?
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param ssc [logical] should the standard errors of the coefficients be corrected for small sample bias?
#' @param df [logical] should the degree of freedoms of the Wald statistic be computed using the Satterthwaite correction?
#' @param ... [logical] argument passed to lower level methods.
#' 
#' @seealso \code{\link{sCorrect}} for more detail about the small sample correction.
#'
#' @details \code{summary2} is the same as \code{summary}
#' except that it first computes the small sample correction (but does not store it).
#' So if \code{summary2} is to be called several times,
#' it is more efficient to pre-compute the quantities for the small sample correction
#' using \code{sCorrect} and then call \code{summary2}.
#'
#' \code{summary2} returns an object with an element \code{table2} containing the estimates, standard errors, degrees of freedom,
#' upper and lower limits of the confidence intervals, test statistics, and p-values.
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
  function(object, ...) UseMethod("summary2")

## * summary2.lm
#' @rdname summary2
#' @method summary2 lm
#' @export
summary2.lm <- function(object, ssc = lava.options()$ssc, df = lava.options()$df,
                        ...){

    object.SSC <- sCorrect(object, ssc = ssc, df = df)
    return(summary2(object.SSC, ...))

}
## * summary2.gls
#' @rdname summary2
#' @method summary2 gls
#' @export
summary2.gls <- summary2.lm

## * summary2.lme
#' @rdname summary2
#' @method summary2 lme
#' @export
summary2.lme <- summary2.lm

## * summary2.lvmfit
#' @rdname summary2
#' @method summary2 lvmfit
#' @export
summary2.lvmfit <- summary2.lm

## * summary2.sCorrect
#' @rdname summary2
#' @export
summary2.sCorrect <- function(object, robust = FALSE, digit = max(3, getOption("digit")),
                              ...){

    ## ** new model parameters
    param <- coef2(object, as.lava = TRUE)
    name.param <- names(param)
    n.param <- length(name.param)

    ## ** new Wald test
    table.all <- compare2(object,
                          linfct = name.param,
                          robust = robust,
                          cluster = NULL,
                          ssc = ssc, df = df,
                          F.test = FALSE,
                          as.lava = FALSE)

    tableS.all <- summary(table.all, test = multcomp::adjusted("none"))$table2
    rownames(tableS.all) <- name.param

    ## ** get summary
    object0 <- object
    class(object0) <- setdiff(class(object0),c("sCorrect"))
    object.summary <- summary(object0, digits = digit, ...)

    ## ** lm
    if(inherits(object,"lm")){
        object.summary$coefficients <- tableS.all[name.param,c("estimate","std.error","statistic","df","p.value"),drop=FALSE]        
        dimnames(object.summary$coefficients) <- list(name.param,
                                                      c("Value","Std.Error","t-value","df","p-value")
                                                      )

        object.summary$residuals <- residuals2(object, type = "response", format = "long")
        object.summary$sigma <- sqrt(tableS.all["sigma2","estimate"])

        object.summary$cov.unscaled <- NULL
        object.summary$fstatistic <- NULL
    }
    ## ** gls / lme
    if(inherits(object,"gls") || inherits(object,"lme")){

        object.summary$tTable <- tableS.all[name.param,c("estimate","std.error","statistic","df","p.value"),drop=FALSE]
        dimnames(object.summary$tTable) <- list(name.param,
                                                c("Value","Std.Error","t-value","df","p-value")
                                                )

        object.summary$residuals <- quantile(residuals2(object, type = "normalized", format = "long"),
                                             na.rm = TRUE)        
        object.summary$sigma <- sqrt(tableS.all["sigma2","estimate"])
    }
    
    ## ** lvmfit
    if(inherits(object,"lvmfit")){
        previous.summary <- object.summary$coef
        object.summary$coef <- tableS.all[name.param,c("estimate","std.error","statistic","df","p.value"),drop=FALSE]

        if(!is.null(object$cluster) || inherits(object,"lvm.missing")){
        
            ## if(robust == FALSE){
            ##     stop("Can only print summary for robust standard errors \n",
            ##          "when the object contain a cluster variable \n")
            ## }
            colnames(object.summary$coef) <- c("Estimate","Std. Error","t-value","df","P-value")
            object.summary$coef[,"t-value"] <- NA

            colnames(object.summary$coefmat) <- c("Estimate","Std. Error","t-value","P-value", "std.xy")
            object.summary$coefmat[,"t-value"] <- ""
        
        }else{
            colnames(object.summary$coef) <- c("Estimate", "Std. Error", "t-value", "df", "P-value")
        }

        ## find digit
        vec.char <- setdiff(object.summary$coefmat[,"Estimate"],"")
        digit <- max(c(nchar(gsub(".","",vec.char,fixed = TRUE)))-1,1)

        ## ** update summary
        ## *** vcov
        object.summary$vcov <- attr(object$dVcov, "vcov.param")[name.param,name.param]    

        ## *** coef
        lava.rownames <- rownames(previous.summary)
        ## add rows corresponding to reference parameters
        missing.rows <- setdiff(lava.rownames,rownames(object.summary$coef))
        if(length(missing.rows)>0){
            addon <- previous.summary[missing.rows,
                                      c("Estimate","Std. Error","Z-value","P-value"),
                                      drop=FALSE]
            colnames(addon)[3] <- "t-value"
            object.summary$coef <- rbind(object.summary$coef, cbind(addon,df=NA))
        }

        ## re-order table according to lava
        object.summary$coef <- object.summary$coef[rownames(previous.summary),,drop=FALSE]
        ## remove unappropriate p.values
        lava.NApvalue <- which(is.na(previous.summary[,"P-value"]))
        object.summary$coef[lava.NApvalue,"P-value"] <- NA
    
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
        table.coef <- object.summary$coef
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
    }

    ## ** gather all results in one table
    object.summary$table2 <- data.frame(matrix(NA, nrow = n.param, ncol = 7,
                                               dimnames = list(name.param,
                                                               c("estimate","std.error","df","ci.lower","ci.upper","statistic","p.value"))
                                               ), stringsAsFactors = FALSE)

    object.summary$table2$estimate <- tableS.all[name.param,"estimate"]
    object.summary$table2$std.error <- tableS.all[name.param,"std.error"]
    object.summary$table2$df <- tableS.all[name.param,"df"]
    object.summary$table2$ci.lower <- object.summary$table2$estimate + object.summary$table2$std.error * qt(p=0.025, df = object.summary$table2$df)
    object.summary$table2$ci.upper <- object.summary$table2$estimate + object.summary$table2$std.error * qt(p=0.975, df = object.summary$table2$df)
    object.summary$table2$statistic <- tableS.all[name.param,"statistic"]
    object.summary$table2$p.value <- tableS.all[name.param,"p.value"]

    ## ** export
    return(object.summary)
    
}

## * confint.sCorrect
confint2 <- confint
confint2.sCorrect <- function(object,...){
    summary2(object)$table2
}

##----------------------------------------------------------------------
### Scorrect-summary2.R ends here


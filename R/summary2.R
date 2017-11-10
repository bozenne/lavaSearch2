### summary2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 10 2017 (10:57) 
## Version: 
## Last-Updated: nov 10 2017 (12:08) 
##           By: Brice Ozenne
##     Update #: 33
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - summary2
#' @title  Summary with small sample correction
#' @description Summary with small sample correction
#' @name summary2
#' @export
`summary2` <-
  function(object, ...) UseMethod("summary2")

## * summary2.lvmfit
#' @rdname summary
#' @export
summary2.lvmfit <- function(object, digits = max(3, getOption("digit")),
                            robust = FALSE,
                            adjust.residuals = FALSE, power = 0.5, as.clubSandwich = TRUE,
                            ...){

    object.summary <- summary(object, digits = digits, ...)
    param <- pars(object)
    name.param <- names(param)
    name.allParam <- rownames(object.summary$coef)
    n.allParam <- length(name.allParam)
    data <- model.frame(object)
    
### ** update variance covariance matrix
    if(robust == FALSE){
        vcov.object <- vcov(object)
        attr(vcov.object, "det") <- NULL
        attr(vcov.object, "pseudo") <- NULL
        attr(vcov.object, "minSV") <- NULL
    }

### ** compute degrees of freedom
    df.adj <- dfVariance(object, p = param, data = data, vcov.param = vcov.object,
                         robust = robust,
                         adjust.residuals = adjust.residuals, power = power, as.clubSandwich = as.clubSandwich, ...)

### ** update summary
    ## vcov
    object.summary$vcov <- vcov.object[name.param,name.param]    

    ## coef
    table.coef <- matrix(NA, ncol = 5, nrow = n.allParam,
                         dimnames = list(name.allParam,
                                         c("Estimate", "Std. Error", "t-value", "df", "P-value")))
    table.coef[,"Estimate"] <- object.summary$coef[,"Estimate"]
    table.coef[name.param,"Std. Error"] <- sqrt(diag(object.summary$vcov))
    table.coef[name.param,"t-value"] <- table.coef[name.param,"Estimate"]/table.coef[name.param,"Std. Error"]
    table.coef[name.param,"df"] <- df.adj
    table.coef[name.param,"P-value"] <- unlist(Map(function(iT, iDF){
        2*(1-pt(abs(iT), df = iDF))
    }, 
    iT = table.coef[name.param,"t-value"],
    iDF = table.coef[name.param,"df"]))

    object.summary$coef <- table.coef

    ## coefmat
    name.label0 <- trimws(rownames(CoefMat(e.lvm, labels = 0, level = 9)), which = "both")
    index.titleVariance <- which(name.label0=="Residual Variances:")
    if(length(index.titleVariance)>0){
        index.titleVariance <- (index.titleVariance+1):length(name.label0)
        name.label0[index.titleVariance] <- paste0(name.label0[index.titleVariance],lava.options()$symbols[2],name.label0[index.titleVariance])
    }

    table.coefmat <- object.summary$coefmat
    colnames(table.coefmat)[3:5] <- c("t-value","df","P-value")

    sd2add <- formatC(table.coef[,"Std. Error"], digit = digits - 1, format = "g",  preserve.width = "common", flag = "")
    sd2add <- gsub(" NA","",sd2add)
    df2add <- as.character(round(table.coef[,"df"],2))    
    df2add[is.na(df2add)] <- ""
    t2add <- formatC(table.coef[,"t-value"], digit = digits - 1, format = "g",  preserve.width = "common", flag = "")
    t2add <- gsub(" NA","",t2add)
    p2add <- formatC(table.coef[,"P-value"], digit = digits - 1, format = "g",  preserve.width = "common", flag = "")
    p2add <- gsub(" NA","",p2add)

    M2add <- cbind(sd2add,df2add,t2add,p2add)
    table.coefmat[match(rownames(table.coef), name.label0),c("Std. Error","t-value","df","P-value")] <- M2add
    object.summary$coefmat <- table.coefmat

    class(object.summary$coefmat[,2])
    
### ** Export
    return(object.summary)    
}

##----------------------------------------------------------------------
### summary2.R ends here

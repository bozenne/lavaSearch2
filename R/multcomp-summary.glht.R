### summary.glht2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  2 2018 (09:20) 
## Version: 
## Last-Updated: jan 24 2020 (15:23) 
##           By: Brice Ozenne
##     Update #: 83
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * summary.glht2
summary.glht2 <- function(object, confint = TRUE, conf.level = 0.95, ...){

    keep.class <- class(object)
    class(object) <- setdiff(keep.class, "glht2")
    output <- summary(object, ...)

    name.hypo <- rownames(output$linfct)
    n.hypo <- length(name.hypo)

    if(confint && is.null(output$confint) && output$test$type %in% c("none","bonferroni","single-step")){
        if(output$test$type == "none"){
            output <- confint(output, level = conf.level, calpha = univariate_calpha())
        }else if(output$test$type == "bonferroni"){
            output <- confint(output, level = 1-(1-conf.level)/n.hypo, calpha = univariate_calpha())
        }else if(output$test$type == "single-step"){
            output <- confint(output, level = conf.level, calpha = adjusted_calpha())
        }else{
            output$confint <- matrix(NA, nrow = n.hypo, ncol = 3,
                                     dimnames = list(name.hypo, c("Estimate","lwr","upr")))
        }
    }

    output$test2 <- data.frame(matrix(NA, nrow = n.hypo, ncol = 7,
                                      dimnames = list(paste0(name.hypo, " == ", output$rhs),
                                                      c("estimate","std.error","df","ci.lower","ci.upper","statistic","p.value"))
                                      ), stringsAsFactors = FALSE)
    output$test2$estimate <- output$test$coefficients
    output$test2$std.error <- output$test$sigma
    output$test2$df <- output$df
    output$test2$df[output$test2$df==0] <- Inf
    output$test2$ci.lower <- output$confint[,"lwr"]
    output$test2$ci.upper <- output$confint[,"upr"]
    output$test2$statistic <- output$test$tstat
    output$test2$p.value <- output$test$pvalues

    class(output) <- append(c("summary.glht2","summary.glht"),keep.class)
    return(output)
}

## * print.summary.glht2
print.summary.glht2 <- function(object, digits = 3, digits.p.value = 4,
                                columns = c("estimate","std.error","df","ci.lower","ci.upper","statistic","p.value"),
                                ...){
    
    columns <- match.arg(columns, choices = c("estimate","std.error","df","ci.lower","ci.upper","statistic","p.value"), several.ok = TRUE)
    type <- object$type
    call <- if(isS4(object$model)){object$model@call}else{object$model$call}
    alternative <- object$alternativ
    type <- object$test$type
    txt.type <- switch(type,
                       "univariate" = "(Univariate p values reported)", 
                       "single-step" = paste0("(Adjusted p values reported -- max-test method (single-step)"), 
                       "free" = paste0("(Adjusted p values reported -- max-test method (step-down)"), 
                       "Westfall" = paste0("(Adjusted p values reported -- max-test method (step-down with logical restrictions))"), 
                       paste0("(Adjusted p values reported --", type, "method)")
                       )
    txt.robust <- switch(as.character(object$robust),
                         "TRUE" = "Robust standard errors",
                         "FALSE" = "Model-based standard errors"
                         )

    txt.correction <- switch(as.character(object$ssc),
                             "Cox" = " corrected for small sample bias (Cox correction)",
                             "residuals" = " corrected for small sample bias (residual correction)",
                             "NA" = ""
                             )
    
    txt.alternative <- switch(alternative,
                              "less" = "one sided tests - inferiority",
                              "greater" = "one sided tests - superiority",
                              "two.sided" = "two sided tests")

    ## display
    cat("\n\t", "Simultaneous Tests for General Linear Hypotheses\n\n")
    if (!is.null(type)) {
        cat("Multiple Comparisons of Means (",txt.alternative,") \n\n", sep = "")
    }
    if (!is.null(call)) {
        cat("Fit: ")
        print(call)
        cat("\n")
    }
    cat("Linear Hypotheses:\n")
    stats::printCoefmat(object$test2[,columns,drop=FALSE], digits = digits,
                        has.Pvalue = "p.value" %in% columns, P.values = "p.value" %in% columns, eps.Pvalue = 10^{-digits.p.value})

    cat(txt.type,"\n")

    error <- attr(object$test$pvalues,"error")
    if(!is.null(error) && error > 1e-12 && "p.value" %in% columns){
        cat("Error when computing the p-value by numerical integration: ",error,"\n",sep="")
    }
    cat("(",txt.robust,")\n",sep="")
    if(nchar(txt.correction)>0){cat("(",txt.correction,")\n",sep="")}
    cat("\n")
    return(invisible(object))
}


######################################################################
### summary.glht2.R ends here

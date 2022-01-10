### sCorrect-summary.glht2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  2 2018 (09:20) 
## Version: 
## Last-Updated: Jan  5 2022 (12:48) 
##           By: Brice Ozenne
##     Update #: 191
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @param transform [function] function to backtransform the estimates and the associated confidence intervals
#' (e.g. \code{exp} if the outcomes have been log-transformed).

## * summary.glht2
#' @export
summary.glht2 <- function(object, confint = TRUE, conf.level = 0.95, transform = NULL, seed = NA, ...){
    if(!is.na(seed)){
        old.seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit( assign(".Random.seed", old.seed, envir = .GlobalEnv, inherits = FALSE) )
        set.seed(seed)
    }
    
    keep.class <- class(object)
    object$test <- NULL
    object$confint <- NULL
    class(object) <- setdiff(keep.class, "glht2")
    keep.df <- object$df
    test.df <- any( (keep.df>0) * (!is.infinite(keep.df)) == 1 )
    object$df <- round(median(object$df))
    output <- summary(object, ...)
    ## restaure df when possible
    method.adjust <- output$test$type
    if(NROW(object$linfct)==1){method.adjust <- "none"}
    if(test.df && method.adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none","univariate")){
        output$df <- keep.df
        output$test$pvalues <- p.adjust(2*(1-pt(abs(output$test$tstat), df = keep.df)), method = method.adjust)
    }
    
    name.hypo <- rownames(output$linfct)
    n.hypo <- length(name.hypo)
    if(confint && method.adjust %in% c("univariate","none","bonferroni","single-step")){
        if(method.adjust %in% c("none","univariate","bonferroni")){
            alpha <- switch(method.adjust,
                            "none" = 1-conf.level,
                            "univariate" = 1-conf.level,
                            "bonferroni" = (1-conf.level)/n.hypo)
            if(test.df){
                q <- qt(1-alpha/2, df = output$df)
            }else{
                q <- qnorm(1-alpha/2)
            }
            output$confint <- data.frame(matrix(NA, ncol = 3, nrow = n.hypo,
                                                dimnames = list(name.hypo, c("Estimate","lwr","upr"))))
            output$confint$Estimate <- as.double(output$test$coef)
            output$confint$lwr <- as.double(output$test$coef - q * output$test$sigma)
            output$confint$upr <- as.double(output$test$coef + q * output$test$sigma)
            ## range(confint(output, level = 1-alpha, calpha = univariate_calpha())$confint-output$confint)
        }else if(method.adjust == "single-step"){
            output <- confint(output, level = conf.level, calpha = adjusted_calpha())
        }else{
            output$confint <- matrix(NA, nrow = n.hypo, ncol = 3,
                                     dimnames = list(name.hypo, c("Estimate","lwr","upr")))
        }
    }
    output$table2 <- data.frame(matrix(NA, nrow = n.hypo, ncol = 7,
                                       dimnames = list(paste0(name.hypo, " == ", output$rhs),
                                                       c("estimate","se","df","lower","upper","statistic","p.value"))
                                       ), stringsAsFactors = FALSE)
    output$table2$estimate <- output$test$coefficients
    output$table2$se <- output$test$sigma
    output$table2$df <- output$df
    output$table2$df[output$table2$df==0] <- Inf
    output$table2$lower <- output$confint[,"lwr"]
    output$table2$upper <- output$confint[,"upr"]
    output$table2$statistic <- output$test$tstat
    output$table2$p.value <- output$test$pvalues
    output$seed <- seed
    
    ## ** transformation
    output$table2 <- transformSummaryTable(output$table2,
                                           transform = transform,
                                           conf.level = conf.level)

    ## ** export    
    class(output) <- append(c("summary.glht2","summary.glht"),keep.class)
    return(output)
}

## * print.summary.glht2
#' @export
print.summary.glht2 <- function(object,
                                digits = max(3L, getOption("digits") - 2L),
                                digits.p.value = max(3L, getOption("digits") - 2L),
                                columns = c("estimate","se","df","lower","upper","statistic","p.value"),
                                ...){
    
    columns <- match.arg(columns, choices = c("estimate","se","df","lower","upper","statistic","p.value"), several.ok = TRUE)
    type <- object$type
    call <- if(isS4(object$model)){object$model@call}else{object$model$call}
    alternative <- object$alternativ
    type <- object$test$type
    txt.type <- switch(type,
                       "univariate" = "(Univariate p values reported)", 
                       "single-step" = paste0("(Adjusted p values reported -- single step max-test)"), 
                       "free" = paste0("(Adjusted p values reported -- step down max-test)"), 
                       "Westfall" = paste0("(Adjusted p values reported -- step down max-test with logical restrictions)"), 
                       paste0("(Adjusted p values reported -- ", type, " method)")
                       )
    txt.robust <- switch(as.character(object$robust),
                         "TRUE" = "Robust",
                         "FALSE" = "Model-based"
                         )

    ## txt.correction <- switch(as.character(object$ssc),
    ##                          "Cox" = " corrected for small sample bias (Cox correction)",
    ##                          "residuals" = " corrected for small sample bias (residual correction)",
    ##                          "NA" = ""
    ##                          )
    
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
        cat("Standard errors: ",txt.robust,"\n",sep="")
        cat("\n")
    }
    cat("Linear Hypotheses:\n")
    stats::printCoefmat(object$table2[,columns[columns %in% names(object$table2)],drop=FALSE], digits = digits,
                        has.Pvalue = "p.value" %in% columns,
                        P.values = "p.value" %in% columns,
                        eps.Pvalue = 10^{-digits.p.value})

    cat(txt.type,"\n")
    error <- attr(object$test$pvalues,"error")
    if(!is.null(error) && error > 1e-12 && "p.value" %in% columns){
        txt.error <- paste0("Error when computing the p-value by numerical integration: ", signif(error, digits = digits))
        if(!is.na(object$seed)){
            txt.error <- paste0(txt.error," (seed ",object$seed,")")
        }
        cat(txt.error,"\n")
    }

    
    if(!is.null(object$global)){
        cat("\nGlobal test: p.value=",format.pval(object$global["p.value"], digits = digits, eps = 10^(-digits.p.value)),
            " (statistic=",round(object$global["statistic"], digits = digits),
            ", df=",round(object$global["df"], digits = digits),")\n",sep="")
    }
    ## if(nchar(txt.correction)>0){cat("(",txt.correction,")\n",sep="")}
    cat("\n")
    return(invisible(object))
}


######################################################################
### sCorrect-summary.glht2.R ends here

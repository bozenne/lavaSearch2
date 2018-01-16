### summary.modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 30 2017 (10:46) 
## Version: 
## last-updated: jan 12 2018 (10:39) 
##           By: Brice Ozenne
##     Update #: 53
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * method summary.modelsearch2
#' @method summary modelsearch2
#' @export
summary.modelsearch2 <- function(object, display = TRUE, ...){

    convergence <- df <- NULL ## [:for CRAN check] data.table
    
    ## ** extract data from object
    xx <- copy(object$sequenceTest)
    dt.seqTest <- rbindlist(lapply(xx, function(step){
        if("convergence" %in% names(step)){
            step[, c("noConvergence") := sum(convergence!=0)]
            step[, c("convergence") := sum(convergence==0)]
        }
        indexMax <- which.max(abs(step$statistic))
        return(step[indexMax])
    }))
    n.step <- NROW(dt.seqTest)
    n.selected <- sum(dt.seqTest$selected)

    keep.cols <- c("link","nTests","noConvergence","statistic","adjusted.p.value")    
    if(!is.na(object$method.p.adjust) && object$method.p.adjust == "max"){
        keep.cols <- c(keep.cols,"quantile")
    }

    if("df" %in% names(dt.seqTest)){
        dt.seqTest[, "nTests.adj" := 0.05/(2*(1-stats::pt(quantile, df = df)))]
    }else{
        dt.seqTest[, "nTests.adj" := 0.05/(2*(1-stats::pnorm(quantile)))]
    }
    
    ## ** output
    out <- list(output = list(), data = dt.seqTest)
    statistic <- switch(object$statistic,
                        "Wald" = "robust Wald",
                        "score" = "score",
                        "LR" = "likelihood ratio",
                        "NA" = "NA")
           
    out$output$message.pre <- paste0("Sequential search for local dependence using the ",statistic," statistic \n")
    if(n.selected==0){
        out$output$message.pre <- c(out$output$message.pre,
                                    "The variable selection procedure did not retain any variable \n")
    }else{
        out$output$message.pre <- c(out$output$message.pre,
                                    paste0("The variable selection procedure retained ",n.selected," variable",
                                           if(n.selected>1){"s"},":\n")
                                    )     
    }
     
    out$output$table <- dt.seqTest[,.SD,.SDcols = keep.cols]
    data.table::setcolorder(out$output$table, keep.cols)
    out$output$message.post <- paste0("confidence level: ",1-object$alpha," (two sided, adjustement: ",object$method.p.adjust,")\n")  

    ## ** display
    if(display){
        cat(out$output$message.pre)
        print(out$output$table)
        cat(out$output$message.post)        
    }
    
    ## ** export
    return(invisible(out))
}



#----------------------------------------------------------------------
### summary.modelsearch2.R ends here

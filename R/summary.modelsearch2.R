### summary.modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 30 2017 (10:46) 
## Version: 
## last-updated: aug 30 2017 (10:56) 
##           By: Brice Ozenne
##     Update #: 14
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @method summary modelsearch2
#' @export
summary.modelsearch2 <- function(object, display = TRUE, ...){

    ## extract data from object
    xx <- copy(object$sequenceTest)
    dt.seqTest <- rbindlist(lapply(xx, function(step){
        if("convergence" %in% names(step)){
            step[, c("convergence") := sum(.SD$convergence==0)]
        }
        indexMax <- which.max(abs(step$statistic))
        return(step[indexMax])
    }))
    n.selected <- sum(dt.seqTest$selected)

    keep.cols <- c("link","statistic","adjusted.p.value","nTests")
    if(object$statistic %in% c("Wald","LR")){
        keep.cols <- c(keep.cols,"convergence")
    }
    if(object$method.p.adjust == "max"){
        keep.cols <- c(keep.cols,"quantile")
    }    

    
    ## output
    out <- list(output = list(), data = dt.seqTest)
    if(n.selected==0){
        out$output$message.pre <- "The variable selection procedure did not retain any variable \n"
    }else{
        out$output$message.pre <- paste0("The variable selection procedure retained ",n.selected," variable",
                                  if(n.selected>1){"s"},":\n")     
    }
    out$output$table <- dt.seqTest[,.SD,.SDcols = keep.cols]

    out$output$message.post <- paste0("confidence level: ",1-object$alpha," (two sided)\n")  

    ## display
    if(display){
        cat(out$output$message.pre)
        print(out$output$table)
        cat(out$output$message.post)        
    }
    
    ## export
    return(invisible(out))
}



#----------------------------------------------------------------------
### summary.modelsearch2.R ends here

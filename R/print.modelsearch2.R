### print.modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 25 2017 (10:18) 
## Version: 
## last-updated: aug 28 2017 (11:40) 
##           By: Brice Ozenne
##     Update #: 9
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ print.modelsearch2
#' @method print modelsearch2
#' @export
print.modelsearch2 <- function(x, ...){

    ## extract data from x
    xx <- copy(x$sequenceTest)
    dt.seqTest <- rbindlist(lapply(xx, function(step){
        if("convergence" %in% names(step)){
            step[, c("convergence") := sum(.SD$convergence==0)]
        }
        indexMax <- which.max(abs(step$statistic))
        return(step[indexMax])
    }))
    n.selected <- sum(dt.seqTest$selected)

    keep.cols <- c("link","statistic","adjusted.p.value","nTests")
    if(x$statistic %in% c("Wald","LR")){
        keep.cols <- c(keep.cols,"convergence")
    }
    if(x$method.p.adjust == "max"){
        keep.cols <- c(keep.cols,"quantile")
    }    
    ## display
    if(n.selected==0){
        cat("The variable selection procedure did not retain any variable \n") 
    }else{
        cat("The variable selection procedure retained ",n.selected," variable",
            if(n.selected>1){"s"},":\n", sep = "")     
    }
    print(dt.seqTest[,.SD,.SDcols = keep.cols])
    
    cat("confidence level: ",1-x$alpha," (two sided)\n",sep="")  
    

}
# }}}

#----------------------------------------------------------------------
### print.modelsearch2.R ends here

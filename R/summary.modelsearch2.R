### summary.modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 30 2017 (10:46) 
## Version: 
## last-updated: sep 19 2018 (12:06) 
##           By: Brice Ozenne
##     Update #: 89
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * method summary.modelsearch2
#' @title summary Method for modelsearch2 Objects
#' @description summary method for modelsearch2 objects.
#'
#' @param object output of the \code{modelsearch2} function.
#' @param print should the summary be printed in the terminal.
#' @param ... [internal] only used by the generic method.
#' 
#' @method summary modelsearch2
#' @export
summary.modelsearch2 <- function(object, print = TRUE, ...){

    p.value <- NULL # [:for CRAN check] subset
    
    ## ** extract data from object
    n.step <- nStep(object)
    tableTest <- do.call(rbind,lapply(object$sequenceTest, function(iTest){
        iTest[which.max(iTest$statistic),]
    }))
    n.selected <- sum(tableTest$selected)

    keep.cols <- c("link","nTests","statistic","adjusted.p.value")

    ## ** output
    out <- list()
    out$message.pre <- paste0("Sequential search for local dependence using the score statistic \n")
    if(n.selected==0){
        out$message.pre <- c(out$message.pre,
                             "The variable selection procedure did not retain any variable \n")
    }else{
        out$message.pre <- c(out$message.pre,
                             paste0("The variable selection procedure retained ",n.selected," variable",
                                    if(n.selected>1){"s"},":\n")
                             )     
    }

    out$table <- tableTest
    out$message.post <- paste0("Confidence level: ",1-object$alpha," (two sided, adjustement: ",object$method.p.adjust,")\n")  

    ## ** display
    if(print){
        cat(out$message.pre,sep="")
        print(out$table)
        cat(out$message.post,sep="")        
    }
    
    ## ** export
    return(invisible(out))
}



#----------------------------------------------------------------------
### summary.modelsearch2.R ends here

### summary.glht2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  2 2018 (09:20) 
## Version: 
## Last-Updated: maj  2 2018 (10:38) 
##           By: Brice Ozenne
##     Update #: 20
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * method summary.glht2
#' @title summary Method for glht2 Objects
#' @description summary method for glht2 objects.
#'
#' @param object output of the \code{glht2} function.
#' @param print should the summary be printed in the terminal.
#' @param ... [internal] only used by the generic method.
#' 
#' @method summary glht2
#' @export
summary.glht2 <- function(object, print = TRUE, ...){

    class(object) <- setdiff(class(object), "glht2")
    object.summary <- summary(object)
    output <- utils::capture.output(print(object.summary))
    
    txt.robust <- switch(as.character(object$robust),
                         "TRUE" = "Robust standard errors",
                         "FALSE" = "Model-based standard errors"
                         )
    txt.correction <- switch(as.character(object$bias.correct),
                             "TRUE" = " corrected for small sample bias",
                             "FALSE" = ""
                             )
    output[length(output)] <- paste0("(",txt.robust,txt.correction,")\n")

    if(print){
        cat(paste0(output,collapse = "\n"))
    }
    
    return(invisible(object.summary))
}


######################################################################
### summary.glht2.R ends here

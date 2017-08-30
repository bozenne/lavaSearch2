### print.modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 25 2017 (10:18) 
## Version: 
## last-updated: aug 30 2017 (10:58) 
##           By: Brice Ozenne
##     Update #: 11
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

    out <- summary(x, display = TRUE, ...)
    return(invisible(out))
}
# }}}

#----------------------------------------------------------------------
### print.modelsearch2.R ends here

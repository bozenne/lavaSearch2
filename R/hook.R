### hook.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 28 2018 (09:53) 
## Version: 
## Last-Updated: nov 28 2018 (09:54) 
##           By: Brice Ozenne
##     Update #: 1
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * hooks
#' @title Hooks for lavaSearch2
#' @description Get and add hook for lavaSeach2
#' @name hook.lavaSearch2
#'
#' @param x the function to add to the hook
#' @param hook the name of the hook
#' @param ... for compatibility with lava

#' @rdname hook.reduce
#' @export
gethook_lavaSearch2<- function (hook, ...){
  get(hook, envir = lavaSearch2.env)
}

#' @rdname hook.reduce
#' @export
addhook_lavaSearch2 <- function (x, hook, ...){
  newhooks <- unique(c(gethook_lavaSearch2(hook), x))
  assign(hook, newhooks, envir = lavaSearch2.env)
  invisible(newhooks)
}



######################################################################
### hook.R ends here

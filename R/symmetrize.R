### symmetrize.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (09:42) 
## Version: 
## Last-Updated: nov  8 2017 (09:49) 
##           By: Brice Ozenne
##     Update #: 8
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @title Symmetrize a matrix
#' @description Complete the upper (or lower) extra-diagonal terms  in order to obtain a symmetric matrix.
#' @param M a matrix
#' @param upper should the upp
#' @examples
#' M <- matrix(NA, 4, 4)
#' M[lower.tri(M)] <- 1:6
#'
#' symmetrize(M, update.upper = TRUE) # good
#'
#' M[upper.tri(M, diag = FALSE)] <- M[lower.tri(M, diag = FALSE)]
#' M # wrong
#' @export
symmetrize <- function(M, update.upper = TRUE){

    if(update.upper){
        M[upper.tri(M, diag = FALSE)] <- t(M)[upper.tri(M, diag = FALSE)]
    }else{
        M[lower.tri(M, diag = FALSE)] <- t(M)[lower.tri(M, diag = FALSE)]
    }
    return(M)    
}

##----------------------------------------------------------------------
### symmetrize.R ends here

### checkData.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 26 2017 (14:25) 
## Version: 
## last-updated: okt 26 2017 (16:34) 
##           By: Brice Ozenne
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - checkData
#' @title Check that validity of the dataset
#' @description Check the validity of the dataset used to estimate a lvm.
#' 
#' @name checkData
#' 
#' @param x a lvm model
#' @param data the dataset containing the variables used to estimate the lvm.
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x
#' latent(m) <- ~u
#'
#' d <- sim(m,1e2)
#'
#' try(checkData(m, d)) # return an error
#' 
#' checkData(m, d[,-4])
#' 
#' try(checkData(m, d[,-(3:4)])) # return an error
#'
#' @export
`checkData` <-
  function(x, data) UseMethod("checkData")

## * checkData.lvm
#' @rdname checkData
#' @export
checkData.lvm <- function(x, data){ 
    vars <- vars(x)
    latent <- latent(x)
    missingVars <- vars[vars %in% names(data) == FALSE]

    if(!identical(sort(latent),sort(missingVars))){
        stop("Wrong specification of the latent variables \n",
             "latent variables according to the LVM: ",paste(latent, collapse = " "),"\n",
             "missing variables in data: ",paste(missingVars, collapse = " "),"\n")
        return(invisible(FALSE))
    }else{
        return(invisible(TRUE))
    }
}


#----------------------------------------------------------------------
### checkData.R ends here

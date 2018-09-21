### methods-modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep 22 2017 (16:43) 
## Version: 
## last-updated: sep 21 2018 (16:43) 
##           By: Brice Ozenne
##     Update #: 159
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * nStep
## ** documentation - nStep
#' @title Find the Number of Steps Performed During the Sequential Testing
#' @description Find the number of steps performed during the sequential testing.
#' @name nStep
#' 
#' @param object a \code{modelsearch2} object.
#' 
#' @return an integer.
#' 
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, method.p.adjust = "holm")
#' nStep(res)
#'
#' @concept modelsearch
#' @concept extractor
#' @export
`nStep` <-
  function(object) UseMethod("nStep")

## ** function - nStep
#' @rdname nStep
#' @export
nStep.modelsearch2 <- function(object){

    return(length(object$sequenceTest))
    
}

## * getStep
## ** documentation - getStep
#' @title Extract one Step From the Sequential Procedure
#' @description Extract one step from the sequential procedure.
#' @name getStep
#' 
#' @param object a \code{modelsearch2} object
#' @param step [integer >0] which test should be extracted?
#' @param slot [character] the element from the modelsearch2 object that should be extracted.
#' @param ... [internal] only used by the generic method.
#'
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, method.p.adjust = "holm")
#' getStep(res)
#' getStep(res, slot = "sequenceTest")
#' getStep(res, slot = "sequenceQuantile")
#' getStep(res, step = 1)
#'
#' @concept modelsearch
#' @concept extractor
#' @export
`getStep` <-
  function(object, ...) UseMethod("getStep")

## ** function - getStep
#' @rdname getStep
#' @export
getStep.modelsearch2 <- function(object, step = nStep(object), slot = NULL, ...){
    
    ## ** normalize arguments
    lastStep <- nStep(object)
    if(any(step %in% 1:lastStep == FALSE)){
        stop("step must be an integer between 1 and ",lastStep,"\n")
    }
    if(!is.null(slot) && slot %in% names(object) == FALSE){
        stop("argument \'slot\' must be one of \"",paste(names(object),collapse="\" \""),"\"\n")
    }

    ## ** subset
    new.object <- object
    new.object$sequenceTest <- object$sequenceTest[step]
    new.object$sequenceModel <- object$sequenceModel[step]

    ## ** export
    if(is.null(slot)){
        return(new.object)
    }else{
        if(is.list(new.object[[slot]])){
            return(new.object[[slot]][[1]])
        }else{
            return(new.object[[slot]])
        }
    }
}

## * getNewLink
## ** documentation - getNewLink
#' @title Extract the Links that Have Been Found by the modelsearch2.
#' @description Extract the links that have been found relevant by modelsearch2.
#' @name getNewLink
#' 
#' @param object a \code{modelsearch2} object.
#' @param step [logical] which test should be extracted?
#' @param ... [internal] only used by the generic method.
#'
#' @return A character vector.
#' 
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' df.data <- lava::sim(mSim, 1e2)
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = df.data)
#' res <- modelsearch2(e.lvm, method.p.adjust = "holm")
#' getNewLink(res)
#'
#' @concept modelsearch
#' @concept extractor
#' 
#' @export
`getNewLink` <-
  function(object, ...) UseMethod("getNewLink")

## ** function - getNewLink
#' @rdname getNewLink
#' @export
getNewLink.modelsearch2 <- function(object, step = 1:nStep(object), ...){

    selected <- link <- NULL
    
    ## ** normalize arguments
    lastStep <- nStep(object)
    if(any(step %in% 1:lastStep == FALSE)){
        stop("step must be an integer between 1 and ",lastStep,"\n")
    }

    ## ** extract
    ls.link <- lapply(step, function(x){
        iStep <- getStep(object,step=x,slot="sequenceTest")
        return(subset(iStep, subset = selected == TRUE, select = "link", drop = TRUE))
    })

    return(unlist(ls.link))    
}

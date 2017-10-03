### methods-modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep 22 2017 (16:43) 
## Version: 
## last-updated: okt  3 2017 (17:54) 
##           By: Brice Ozenne
##     Update #: 85
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
#' @title Find the number of steps performed during the sequential testing
#' @description Find the number of steps performed during the sequential testing
#' @name nStep
#' 
#' @param object a modelsearch2 object
#' @examples
#'  mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' dt <- as.data.table(sim(mSim, 1e2))
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = dt)
#' res <- modelsearch2(e.lvm, statistic = "score", method.p.adjust = "holm")
#' nStep(res)
#'
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
#' @title Extract one step from the sequential procedure
#' @description Extract one step from the sequential procedure
#' @name getStep
#' 
#' @param object a modelsearch2 object
#' @param step which test should be extracted?
#' @param slot the element from the modelsearch2 object that should be extracted.
#' @param ... not used.
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' dt <- as.data.table(sim(mSim, 1e2))
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = dt)
#' res <- modelsearch2(e.lvm, statistic = "score", method.p.adjust = "holm")
#' getStep(res)
#' getStep(res, slot = "sequenceTest")
#' getStep(res, slot = "sequenceQuantile")
#' getStep(res, step = 1)
#'
#' @export
`getStep` <-
  function(object, ...) UseMethod("getStep")

## ** function - getStep
#' @rdname getStep
#' @export
getStep.modelsearch2 <- function(object, step = nStep(object), slot = NULL, ...){

    ## ** normalize arguments
    lastStep <- nStep(object)
    if(step %in% 1:lastStep == FALSE){
        stop("step must be an integer between 1 and ",lastStep,"\n")
    }
    if(!is.null(slot) && slot %in% names(object) == FALSE){
        stop("argument \'slot\' must be one of \"",paste(names(object),collapse="\" \""),"\"\n")
    }

    ## ** subset
    new.object <- object
    new.object$sequenceTest <- object$sequenceTest[step]
    new.object$sequenceModel <- object$sequenceModel[step]
    if(object$method.p.adjust == "max"){
        new.object$sequenceQuantile <- object$sequenceQuantile[[step]]
        sequenceIID <- object$sequenceIID[step]
        sequenceSigma <- object$sequenceSigma[step]
    }

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
#' @title Find the links that should be added accroding to the sequential testing
#' @description Find the links that should be added accroding to the sequential testing
#' @name getNewLink
#' 
#' @param object a modelsearch2 object
#' @param step which test should be extracted?
#' @param ... not used
#' 
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' dt <- as.data.table(sim(mSim, 1e2))
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = dt)
#' res <- modelsearch2(e.lvm, statistic = "score", method.p.adjust = "holm")
#' getNewLink(res)
#'
#' @export
`getNewLink` <-
  function(object, ...) UseMethod("getNewLink")

## ** function - getNewLink
#' @rdname getNewLink
#' @export
getNewLink.modelsearch2 <- function(object, step = nStep(object), ...){

    selected <- link <- NULL
    
    ## ** normalize arguments
    lastStep <- nStep(object)
    if(step %in% 1:lastStep == FALSE){
        stop("step must be an integer between 1 and ",lastStep,"\n")
    }

    ## ** extract
    ls.link <- lapply(1:step, function(x){
        getStep(object,step=x,slot="sequenceTest")[selected == TRUE,link]
    })

    return(unlist(ls.link))    
}

## * merge
## ** documentation - merge
#' @title Merge two modelsearch objects
#' @description Merge two modelsearch objects. Does not check for meaningful result.
#' @name merge
#' 
#' @param x a modelsearch2 object.
#' @param y a modelsearch2 object that will be added to x.
#' @param ... not used.
#' 
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' dt <- as.data.table(sim(mSim, 1e2))
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = dt)
#' res.x <- modelsearch2(e.lvm, statistic = "score", method.p.adjust = "holm", nStep = 2)
#' res.y <- modelsearch2(getStep(res.x, slot = "sequenceModel"), statistic = "score", method.p.adjust = "holm")
#' res.xy <- merge(res.x,res.y)
#'
#' modelsearch2(e.lvm, statistic = "score", method.p.adjust = "holm")

## ** function - merge
#' @rdname merge
#' @export
merge.modelsearch2 <- function(x, y, ...){

    ## ** merge
    x$sequenceTest <- c(x$sequenceTest,y$sequenceTest)
    x$sequenceModel <- c(x$sequenceTest,y$sequenceModel)
                
    if(sum(c(y$method.p.adjust,x$method.p.adjust) == "max") == 1){
        if(x$method.p.adjust != "max"){
            lastStep.x <- nStep(x)
            Mtest.x <- getStep(x, step = 1, slot = "sequenceTest")
            nLink.x <- NROW(Mtest.x)
            name.x <- Mtest.x[["link"]]

            x$sequenceQuantile <- matrix(NA, nrow = nLink.x, ncol = lastStep.x,
                                         dimnames = list(name.x, rep("step",lastStep.x)))            
            x$sequenceIID <- vector(mode = "list", length = lastStep.x)
            x$sequenceSigma <- vector(mode = "list", length = lastStep.x)            
        }
        if(y$method.p.adjust != "max"){
            lastStep.y <- nStep(y)
            Mtest.y <- getStep(y, step = 1, slot = "sequenceTest")
            nLink.y <- NROW(Mtest.y)
            name.y <- Mtest.y[["link"]]

            y$sequenceQuantile <- matrix(NA, nrow = nLink.y, ncol = lastStep.y,
                                         dimnames = list(name.y, rep("step",lastStep.y)))            
            y$sequenceIID <- vector(mode = "list", length = lastStep.y)
            y$sequenceSigma <- vector(mode = "list", length = lastStep.y)            
        }
        
        if(NROW(y$sequenceQuantile)>NROW(x$sequenceQuantile)){
            Mall <- matrix(NA, nrow = NROW(y$sequenceQuantile), ncol = NCOL(x$sequenceQuantile),
                           dimnames = list(rownames(y$sequenceQuantile), colnames(x$sequenceQuantile)))
            Mall[rownames(x$sequenceQuantile),] <- x$sequenceQuantile
            x$sequenceQuantile <- Mall
        }else if(NROW(x$sequenceQuantile)>NROW(y$sequenceQuantile)){
            Mall <- matrix(NA, nrow = NROW(x$sequenceQuantile), ncol = NCOL(y$sequenceQuantile),
                           dimnames = list(rownames(x$sequenceQuantile), colnames(y$sequenceQuantile)))
            Mall[rownames(y$sequenceQuantile),] <- y$sequenceQuantile
            y$sequenceQuantile <- Mall            
        }
        
        x$sequenceQuantile <- cbind(x$sequenceQuantile,y$sequenceQuantile)
        x$sequenceIID <- c(x$sequenceIID,y$sequenceIID)
        x$sequenceSigma <- c(x$sequenceSigma,y$sequenceSigma)

    }
    
    for(iSlot in c("statistic","method.p.adjust","alpha","method.iid","cv")){
        x[[iSlot]] <- unique(c(x[[iSlot]],y[[iSlot]]))
        if(length(x[[iSlot]])>1){x[[iSlot]] <- NA}
    }
    ## ** export
    return(x)    
}


#----------------------------------------------------------------------
### methods-modelsearch2.R ends here

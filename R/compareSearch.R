### compareSearch.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: sep 22 2017 (11:57) 
## Version: 
## last-updated: okt  5 2017 (12:42) 
##           By: Brice Ozenne
##     Update #: 166
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - compareSearch
#' @title Compare methods to identify missing local dependencies in a LVM
#' @description Compare methods to identify missing local dependencies in a LVM
#' @name compareSearch
#' 
#' @param object a lvm model
#' @inherit modelsearch2
#' 
#' @examples
#' mSim <- lvm(Y~G+X1+X2)
#' addvar(mSim) <- ~Z1+Z2+Z3+Z4+Z5+Z6
#' dt <- as.data.table(lava::sim(mSim, 1e2))
#'
#' mBase <- lvm(Y~G)
#' addvar(mBase) <- ~X1+X2+Z1+Z2+Z3+Z4+Z5+Z6
#' e.lvm <- estimate(mBase, data = dt)
#'
#' res <- modelsearch2(e.lvm, statistic = "Wald", method.p.adjust = "none")
#' res$sequenceTest
#' \dontrun{
#' res <- compareSearch(e.lvm, statistic = c("score","Wald"),
#'                      method.iid = "iid",
#'                      method.p.adjust = c("holm","fdr","max"))
#' }
#' res
#' 

## * function - compareSearch
#' @rdname compareSearch
#' @export
compareSearch <- function(object, alpha = 0.05,
                          method.iid, method.p.adjust, statistic,
                          trace = 1, ...){

    p.value <- link <- NULL
    
### ** normalize arguments
    method.p.adjust <- sapply(method.p.adjust, match.arg, choices = lava.options()$search.p.adjust, several.ok = TRUE)
    statistic <-  sapply(statistic, match.arg, choices = lava.options()$search.statistic, several.ok = TRUE)
    if("Wald" %in% method.p.adjust){
        method.iid <- match.arg(method.iid, lava.options()$search.iid, several.ok = FALSE)
    }

### ** fit all models using unadjusted p.values
    ls.search <- list()
    if("score" %in% statistic){
        if(trace){
            cat("modelsearch with the score statistic")
        }
        ls.search$score <- modelsearch2(object, statistic = "score", method.p.adjust = "none", method.iid = "iid",
                                        trace = trace-1, ...)
        if(trace){
            cat(" - done \n")
        }
    }
    if("LR" %in% statistic){
        if(trace){
            cat("modelsearch with the log likelihood ratio statistic")
        }
        ls.search$LR <- modelsearch2(object, statistic = "LR", method.p.adjust = "none", method.iid = "iid",
                                     trace = trace-1, ...)
        if(trace){
            cat(" - done \n")
        }
    }
    if("Wald" %in% statistic){
        if(trace){
            cat("modelsearch with the robust Wald statistic")
        }
        if("max" %in% method.p.adjust){
            ls.search$Wald <- modelsearch2(object, statistic = "Wald", method.p.adjust = "max", method.iid = method.iid,
                                           trace = trace-1, ...)

            currentStep <- nStep(ls.search$Wald)
            vec.tempo <- getStep(ls.search$Wald, step = currentStep, slot = "sequenceTest")
            maxStep <- list(...)$nStep
            if(is.null(maxStep)){maxStep <- Inf}
            if(any(vec.tempo$p.value < alpha) && (vec.tempo$selected==FALSE) && (currentStep<maxStep) ){ # continue the modelsearch

                ## ** add the link of the last test to the model (avoid to repeat step)
                model.tempo <- getStep(ls.search$Wald, step=nStep(ls.search$Wald), slot = "sequenceModel")
                link.tempo <- getStep(ls.search$Wald, step=nStep(ls.search$Wald), slot = "sequenceTest")
                newLink.tempo <- link.tempo[which.min(p.value),link]

                ls.args <- lapply(model.tempo$call[-(1:2)], evalInParentEnv)
                restricted.tempo <- unlist(initVarLink(newLink.tempo))
                directive.tempo <- length(grep(lava.options()$symbols[2],newLink.tempo,fixed=TRUE))==0
                
                if("lvmfit" %in% class(object)){
                    model.tempo2 <- .updateModelLink.lvm(model.tempo, args = ls.args,
                                                         restricted = restricted.tempo,
                                                         directive = directive.tempo)
                }else{
                    model.tempo2 <- .updateModelLink.default(model.tempo, args = ls.args,
                                                             restricted = restricted.tempo,
                                                             directive = directive.tempo)
                }
                
                dots <- list(...)
                dots$nStep <- maxStep-currentStep
                if("link" %in% names(dots)){
                    dots$link <- setdiff(dots$link)
                }
                otherSearch <- modelsearch2(model.tempo2, statistic = "Wald", method.p.adjust = "none", method.iid = method.iid,
                                            trace = trace-1, dots)
                
                ls.search$Wald <- merge(ls.search$Wald, otherSearch)   
            }
               
            
        }else{            
            ls.search$Wald <- modelsearch2(object, statistic = "Wald", method.p.adjust = "none", method.iid = method.iid,
                                           trace = FALSE, ...)
        }
        if(trace){
            cat(" - done \n")
        }
    }
    

### ** Adjust p.values
    ls.searchAll <- list()    
    for(iStatistic in statistic){ # iStatistic <- statistic[1]
        for(iAdjust in method.p.adjust){ # iAdjust <- method.p.adjust[1]
            if(iAdjust == "max" && iStatistic != "Wald"){next}
            list.tempo <- list(.adjustModelSearch(ls.search[[iStatistic]], method.p.adjust = iAdjust, alpha  = alpha))
            names(list.tempo) <- paste0(iStatistic,"-",iAdjust)
            ls.searchAll <- c(ls.searchAll,list.tempo)
        }
            
    }
   
### ** Merge results
    ## newlinks
    ls.newlinks <- lapply(ls.searchAll,getNewLink)

    ## value of all links
    name.alllinks <- unique(unlist(lapply(ls.searchAll, function(x){
        names(coef(getStep(x, step = nStep(x), slot = "sequenceModel")))
    })))
    name.search <- names(ls.searchAll)
    table.alllinks <- matrix(NA, nrow = length(name.alllinks), ncol = length(ls.searchAll)+1,
                             dimnames = list(name.alllinks, c("base",name.search)))
    table.alllinks[names(coef(object)),"base"] <- coef(object)
    for(iSearch in name.search){ # iSearch <- name.search[2]
        M.tempo <- getStep(ls.searchAll[[iSearch]], step = nStep(ls.searchAll[[iSearch]]), slot = "sequenceModel")
        table.alllinks[names(coef(M.tempo)),iSearch] <- coef(M.tempo)
    }
    
### ** export
    return(list(newlinks = ls.newlinks,
                table.coef = table.alllinks,
                ls.search = ls.searchAll))
    
}

## * adjustModelSearch
.adjustModelSearch <- function(object, method.p.adjust, alpha){

    object <- copy(object)
    
    ## ** adjust p.value
    seqP.value <- sapply(object$sequenceTest, function(x){        
        if(method.p.adjust!="max"){            
            x[,c("adjusted.p.value") := p.adjust(x$p.value, method = method.p.adjust)]
        }
        return(min(x$adjusted.p.value))
    })

    ## ** stop search when necessary
    index.keepTest <- union(1, which(seqP.value<alpha)+1)
    index.keepTest <- index.keepTest[index.keepTest<=nStep(object)] # remove extra step due to early stop
    
    object$sequenceTest <- object$sequenceTest[index.keepTest]
    object$sequenceModel <- object$sequenceModel[index.keepTest]
    if(method.p.adjust=="max"){ # max never activated
        object$sequenceQuantile <- object$sequenceQuantile[index.keepTest]
        object$sequenceIID <- object$sequenceIID[index.keepTest]
        object$sequenceSigma <- object$sequenceSigma[index.keepTest]
    }

    ## ** update final model
    index.finalModel <- tail(which(seqP.value<alpha),1)
    object$sequenceModel[[length(object$sequenceModel)]] <- object$sequenceModel[[index.finalModel]]

    ## ** update adjustement
    object$method.p.adjust <- method.p.adjust
    
    ## ** export
    return(object)
}

#----------------------------------------------------------------------
### compareSearch.R ends here

### defineCategoricalLink.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 26 2017 (16:35) 
## Version: 
## last-updated: jan 15 2018 (11:41) 
##           By: Brice Ozenne
##     Update #: 119
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - defineCategoricalLink
#' @title Identify Categorical Links in LVM
#' @description Identify categorical links in latent variable models.
#' @name defineCategoricalLink
#' 
#' @param object a lvm model.
#' @param link the links to be analyzed. If NULL, all the coefficients from the lvm model are used instead.
#' @param data the dataset that will be used to fit the model. If \code{NULL}, a simulated data will be generated from the model.
#' 
#' @examples  
#' m <- lvm(Y1~X1+X2,Y2~X1+X3)
#' categorical(m, K = 3) <- "X1"
#' try(defineCategoricalLink(m)) # error
#'
#' categorical(m, K = 3, labels = 1:3) <- "X1"
#' defineCategoricalLink(m)
#' defineCategoricalLink(m, "Y~X1")
#' defineCategoricalLink(m, "X1:0|1")
#' defineCategoricalLink(m, "X1:1|2")
#' defineCategoricalLink(m, c("X1:0|1", "X1:1|2"))
#' defineCategoricalLink(m, c("Y~X1","Y~X2"))
#' defineCategoricalLink(m, c("Y~X2","Y~X1"))
#'
#' @export
`defineCategoricalLink` <-
  function(object, link, data) UseMethod("defineCategoricalLink")


## * defineCategoricalLink.lvm
#' @rdname defineCategoricalLink
#' @export
defineCategoricalLink.lvm <- function(object, link = NULL, data = NULL){

    ### ** normalize arguments
    if(is.null(link)){
        link <- stats::coef(object)
    }
    if(is.null(data)){
        data <- sim(object, 1)
    }
    
    ### ** identify the type of regression variable (continuous or categorical)
    index.cat <- which(link %in% unlist(object$attributes$ordinalparname))
    index.Ncat <- setdiff(1:length(link), index.cat)
    link.Ncat <- setdiff(link[index.Ncat], names(object$attributes$ordinalparname))

    ### ** caracterize links involving categorical variables    
    if(length(index.cat)>0){
        link.cat <- link[index.cat]
        xCAT <- lava_categorical2dummy(object, data)$x

        ## *** find exogenous variable
        X.name.allcat <- sapply(link.cat, function(iL){
            test <- unlist(lapply(object$attributes$ordinalparname, function(vec.coef){iL %in% vec.coef}))
            return(names(object$attributes$ordinalparname)[test])
        })
        UX.name.cat <- unique(X.name.allcat)
    
        ## *** find the level of the exogenous variable
        X.level.cat <- unlist(lapply(UX.name.cat, function(iL){ 
            if(iL %in% names(xCAT$attributes$labels)){
                labels <- xCAT$attributes$labels[[iL]]
                index.label <- which(object$attributes$ordinalparname[[iL]] %in% link.cat)                
                return(labels[1+index.label])
            }else {
                stop("Categorical variables must have labels. Specify argument \'labels\' when calling categorical. \n")
            }            
        }))

        ## *** find endogenous variable
        M.link <- xCAT$M[paste0(X.name.allcat,X.level.cat),,drop = FALSE]
        M.link <- cbind(M.link, as.numeric(rowSums(M.link)==0))
        indexLink <- which(M.link==1, arr.ind = TRUE)
        Y.name.allcat <- colnames(M.link)[indexLink[,"col"]]
            
        ## *** characterize all links
        Xcat.name.allcat <- rownames(M.link)[indexLink[,"row"]]
        X.level.allcat <- as.character(X.level.cat[indexLink[,"row"]])
        external.link.allcat <- link[index.cat[indexLink[,"row"]]]
        original.link.allcat <- paste0(Y.name.allcat, lava.options()$symbol[1], X.name.allcat)
        original.link.allcat[Y.name.allcat == ""] <- gsub("~","",original.link.allcat[Y.name.allcat == ""])
        cat.link.allcat <- paste0(Y.name.allcat, lava.options()$symbol[1], Xcat.name.allcat)
        cat.link.allcat[Y.name.allcat == ""] <- gsub("~","",cat.link.allcat[Y.name.allcat == ""])
        
        dt.cat <- data.table::data.table(link = cat.link.allcat,
                                         endogenous = Y.name.allcat,
                                         exogenous = X.name.allcat,
                                         type = "categorical",
                                         factice = FALSE,
                                         level = X.level.allcat,
                                         originalLink = original.link.allcat,
                                         externalLink = external.link.allcat)       

    }else{
            dt.cat <- NULL
    }

### ** caracterize links involving continuous variables    
    if(length(index.Ncat)>0){

        var.tempo <- initVarLinks(link.Ncat)
        Y.name.Ncat <- var.tempo$var1
        X.name.Ncat <- var.tempo$var2
        test.factice <- X.name.Ncat %in% names(object$attributes$ordinalparname)

        dt.Ncat <- data.table::data.table(link = link.Ncat,
                                          endogenous = Y.name.Ncat,
                                          exogenous = X.name.Ncat,
                                          type = "continuous",
                                          factice = test.factice,
                                          level = NA,
                                          originalLink = link.Ncat,
                                          externalLink = NA)       
    }else{
        dt.Ncat <- NULL
    }

    
    ### ** export
    out <- rbind(dt.Ncat,dt.cat)
    return(out)
}

#----------------------------------------------------------------------
### defineCategoricalLink.R ends here

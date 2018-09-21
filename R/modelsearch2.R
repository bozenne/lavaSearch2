## * modelsearch2 (documentation)
#' @title Data-driven Extension of a Latent Variable Model
#' @description Procedure adding relationship between variables that are supported by the data.
#' @name modelsearch2
#' 
#' @param object a \code{lvmfit} object.
#' @param data [data.frame, optional] the dataset used to identify the model
#' @param link [character, optional for \code{lvmfit} objects] the name of the additional relationships to consider when expanding the model. Should be a vector containing strings like "Y~X". See the details section.
#' @param method.p.adjust [character] the method used to adjust the p.values for multiple comparisons.
#' Can be any method that is valid for the \code{stats::p.adjust} function (e.g. \code{"fdr"}), or \code{"max"} or \code{"fastmax"}.
#' @param alpha [numeric 0-1] the significance cutoff for the p-values.
#' When the p-value is below, the corresponding link will be added to the model
#' and the search will continue. Otherwise the search will stop.
#' @param nStep the maximum number of links that can be added to the model.
#' @param na.omit should tests leading to NA for the test statistic be ignored. Otherwise this will stop the selection process.
#' @param trace [logical] should the execution of the function be traced?
#' @param cpus the number of cpus that can be used for the computations.
#' @param ... additional arguments to be passed to \code{\link{findNewLink}} see details.
#'
#' @details
#' method.p.adjust = \code{"max"} computes the p-values based on the distribution of the max statistic.
#' This max statistic is the max of the square root of the score statistic.
#' The p-value are computed integrating the multivariate normal distribution.
#' 
#' method.p.adjust = \code{"fastmax"} only compute the p-value for the largest statistic.
#' It is faster than \code{"max"} and lead to identical results.
#' 
#' @return A list containing:
#' \itemize{
#' \item sequenceTest: the sequence of test that has been performed.
#' \item sequenceModel: the sequence of models that has been obtained.
#' \item sequenceQuantile: the sequence of rejection threshold. Optional. 
#' \item sequenceIID: the influence functions relative to each test. Optional. 
#' \item sequenceSigma: the covariance matrix relative to each test. Optional. 
#' \item statistic: the argument \code{statistic}.
#' \item method.p.adjust: the argument \code{method.p.adjust}.
#' \item alpha: [numeric 0-1] the significance cutoff for the p-values.
#' \item cv: whether the procedure has converged.
#' } 
#' 
#' @examples
#'
#' #### LVM ####
#' mSim <- lvm()
#' regression(mSim) <- c(y1,y2,y3)~u
#' regression(mSim) <- u~x1+x2
#' categorical(mSim,labels=c("A","B","C")) <- "x2"
#' latent(mSim) <- ~u
#' covariance(mSim) <- y1~y2
#' transform(mSim, Id~u) <- function(x){1:NROW(x)}
#' df.data <- lava::sim(mSim, n = 1e2, latent = FALSE)
#' 
#' m <- lvm(c(y1,y2,y3)~u)
#' latent(m) <- ~u
#' addvar(m) <- ~x1+x2 
#'
#' e <- estimate(m, df.data)
#'
#' resSearch <- modelsearch(e)
#' resSearch
#' resSearch2 <- modelsearch2(e)
#' resSearch2
#'
#' @concept modelsearch
#' @export
`modelsearch2` <-
  function(object, ...) UseMethod("modelsearch2")

## * modelsearch2.lvmfit (code)
#' @rdname modelsearch2
#' @export
modelsearch2.lvmfit <- function(object, link = NULL, data = NULL, 
                                method.p.adjust = "fastmax", alpha = 0.05, 
                                nStep = NULL, na.omit = TRUE, 
                                trace = TRUE, cpus = 1,  
                                ...){

    ## ** check arguments
    ## methods
    method.p.adjust <- match.arg(method.p.adjust, lava.options()$search.p.adjust)    

    ## cpus
    if(is.null(cpus)){ cpus <- parallel::detectCores()}

        if(is.null(cpus) || cpus > 1){
        test.package <- try(requireNamespace("foreach"), silent = TRUE)
        if(inherits(test.package,"try-error")){
            stop("There is no package \'foreach\' \n",
                 "This package is necessary when argument \'cpus\' is greater than 1 \n")
        }
    }
    
    if(!is.null(cpus) && cpus>1){
        if(cpus > parallel::detectCores()){
            stop("Argument \'cpus\' is greater than the number of available CPU cores \n",
                 "available CPU cores: ",parallel::detectCores(),"\n")
        }
    }

    ## extra arguments 
    dots <- list(...)
    if(length(dots)>0){
        stop("modelsearch2 does not take any extra arguments \n",
             "name of the extra arguments: \"",paste(names(dots), collapse = "\" \""),"\" \n")
    }

    ## ** prepare
    ## *** data
    if(is.null(data)){
        data <- as.data.frame(stats::model.frame(object, all = TRUE))
    }
    
    ## *** normalize the links
    if(is.null(link)){
        res.find <- do.call(findNewLink,
                            args = c(list(object$model,
                                          data = data,
                                          output = "names")))
        directive <- res.find$directional
        restricted <- res.find$M.links
        link <- res.find$link
        if(is.null(link)){
            stop("Automatic search has not found any possible additional link \n",
                 "Consider specifying manually the argument \'link\' \n")
        }
    }else{
        resLink <- .initializeLinks(object, data = data, link = link)
        object <- resLink$object
        link <- resLink$link
        directive <- resLink$directive
        restricted <- resLink$restricted
    }    

    ## ** initialization
    if(is.null(nStep)){
        nStep <- NROW(restricted)
    }
    iStep <- 1
    iRestricted <- restricted
    iDirective <- directive
    iLink <- link
    iObject <- object

    ## update of the model
    add.args <- setdiff(names(object$call), c("","object","data","control"))
    ls.call <- lapply(add.args, function(arg){object$call[[arg]]})
    names(ls.call) <- add.args

    ls.call$data <- data
    if(!is.null(data)){
        index.cols <- which(names(data)%in%names(ls.call$data)==FALSE)
        if(length(index.cols)>0){
            ls.call$data <- cbind(ls.call$data,
                                  subset(as.data.frame(data), select = index.cols))
        }
    }
    ls.call$control <- object$control
    ls.call$control$trace <- FALSE
    
    ## output
    ls.seqTests <- list()
    ls.seqModels <- list()
    ls.seqIID <- list() # only for method.p.adjust = "max"
    ls.seqSigma <- list() # only for method.p.adjust = "max"
    vec.seqQuantile <- NULL # only for method.p.adjust = "max"
    
    ## criterion    
    cv <- FALSE

    ## define cluster
    if(cpus>1){
        if(trace>0){
            cl <- parallel::makeCluster(cpus, outfile = "")
        }else{
            cl <- parallel::makeCluster(cpus)
        }
        doParallel::registerDoParallel(cl)

        vec.packages <- c("lavaSearch2")
        parallel::clusterCall(cl, fun = function(x){
            sapply(vec.packages, function(iP){
                suppressPackageStartupMessages(attachNamespace(iP)) ## requireNamespace did not worked
            })
        })
        
    }else{
        cl <- NULL
    }

    ## ** display a summary of the call
    if(trace>0){

        cat("\n",
            "** Sequential variable selection using the score statistic ** \n",
            " Number of possible additional links         : ",length(link)," \n",
            " Maximum number of steps                     : ",nStep,"\n",        
            " Adjustment method for multiple comparisons  : ",method.p.adjust,"\n",
            " Confidence level                            : ",1-alpha,"\n",
            " Number of cpus                              : ",cpus,"\n\n",
            sep="")
    }

    ## ** Forward search
    while(iStep <= nStep && NROW(iRestricted)>0 && cv==FALSE){
        if(trace >= 1){cat("Step ",iStep,":\n",sep="")}

        resStep <- .oneStep_scoresearch(iObject, data = data,
                                        restricted = iRestricted, link = iLink, directive = iDirective,
                                        method.p.adjust = method.p.adjust, 
                                        cl = cl, trace = trace)
        
        ## ** update according the most significant p.value
        ## *** check convergence
        if(na.omit || method.p.adjust == "fastmax"){
            cv <- all(stats::na.omit(resStep$test$adjusted.p.value) > alpha)
            test.na <- FALSE
        }else{
            cv <- all(resStep$test$adjusted.p.value > alpha)
            if(is.na(cv)){                    
                cv <- TRUE
                test.na <- TRUE
            }else{
                test.na <- FALSE
            }
        }

        ## *** identify most promising test
        index.maxTest <- which.max(abs(resStep$test$statistic))[1]
        resStep$test$selected <- FALSE
        resStep$test[index.maxTest,"selected"] <- (resStep$test[index.maxTest,"adjusted.p.value"] <= alpha)
        resStep$test$nTests <- NROW(resStep$test)
        resStep$test <- resStep$test[order(resStep$test$statistic),]

        ## *** update the model
        if(cv==FALSE){
            ls.call$x <- addLink(iObject$model, var1 = iRestricted[index.maxTest,1], var2 = iRestricted[index.maxTest,2],
                                 covariance = 1-iDirective[index.maxTest])

            ## first attempt
            ls.call$start <- stats::coef(iObject)
            suppressWarnings(
                iObject <- tryCatch(do.call(lava::estimate, args = ls.call),
                                      error = function(x){x},
                                      finally = function(x){x})
            )
        
            ## second attempt
            if(inherits(iObject,"error") || iObject$opt$convergence>0){
                ls.call$control$start <- NULL
                suppressWarnings(
                    iObject <- do.call(lava::estimate, args = ls.call)
                )
            }

            ## update links
            iLink <- iLink[-index.maxTest]
            iRestricted <- iRestricted[-index.maxTest,,drop=FALSE]
            iDirective <- iDirective[-index.maxTest]            
        }
        


        ## *** update the output
        ls.seqTests[[iStep]] <- resStep$test
        ls.seqModels[[iStep]] <- iObject
        if(method.p.adjust %in% c("max","fastmax")){
            ls.seqIID[[iStep]] <- resStep$iid
            ls.seqSigma[[iStep]] <- resStep$Sigma
            vec.seqQuantile <- c(vec.seqQuantile,resStep$test$quantile[1])
        }


        ## *** display results
        if(trace > 0){
            rowSelected <- NROW(resStep$test)
            if(cv==FALSE){
                cat("add ",as.character(resStep$test[rowSelected, "link"]),
                    " (statistic = ",resStep$test[rowSelected, "statistic"],
                    ", adjusted.p.value = ",resStep$test[rowSelected, "adjusted.p.value"],
                    ")\n",sep="")
            }else{
                if(test.na){
                    cat("NA among the test statistics \n")
                }else{
                    cat("no variable to add",
                        " (statistic = ",resStep$test[rowSelected, "statistic"],
                        ", adjusted.p.value = ",resStep$test[rowSelected, "adjusted.p.value"],
                        ")\n",sep="")
                }
            }
        }
        iStep <- iStep + 1
    }

    if(cpus>1){
        parallel::stopCluster(cl)
    }

    ## ** export
    if(length(ls.seqIID)==0){ls.seqIID <- NULL}
    if(length(ls.seqSigma)==0){ls.seqSigma <- NULL}
    output <- list(sequenceTest = ls.seqTests,
                   sequenceModel = ls.seqModels,
                   sequenceQuantile = vec.seqQuantile,
                   sequenceIID = ls.seqIID,
                   sequenceSigma = ls.seqSigma,
                   method.p.adjust = method.p.adjust,
                   alpha = alpha,
                   cv = cv)
    class(output) <- "modelsearch2"
    return(output)
}

## * .initializeLinks
.initializeLinks <- function(object, data, link){
    restricted <- do.call(cbind,initVarLinks(link))
    directive <- rep(TRUE, length(link))
    index.Ndir <- grep(lava.options()$symbols[2],link,fixed=TRUE)
    if(length(index.Ndir)>0){
        directive[index.Ndir] <- FALSE
    }

    ## ** get all vars
    if(is.null(data)){            
        data <- evalInParentEnv(object$call$data, envir = environment())             
        if(is.null(data)){
            data <- lava::sim(object,1)
        }            
    }
    ## ** take care of categorical variables
    ls.linkvar <- do.call(rbind,lapply(1:NROW(restricted), function(row){
        data.frame(Y = restricted[row,1],
                   X = var2dummy(object$model,
                                 data = data,
                                 var = restricted[row,2]),
                   dir = directive[row],
                   stringsAsFactors = FALSE)
    }))
    restricted2 <- as.matrix(ls.linkvar[,1:2,drop=FALSE])
    directive <- ls.linkvar[,3]
    link <- paste0(restricted2[,1],lava.options()$symbols[2-directive],restricted2[,2])

    ## ** check links
    allVars.link <- unique(as.vector(restricted2))
    allVars.model <- vars(object$model)
    allVars.data <- names(data)
        
    if(any(allVars.link %in% allVars.model == FALSE)){
        missing.var <- allVars.link[allVars.link %in% allVars.model == FALSE]
            
            if(any(allVars.link %in% allVars.data == FALSE)){
                missing.var <- allVars.link[allVars.link %in% allVars.data == FALSE]
                stop("Some links contains variables that are not in the latent variable model \n",
                     "variables(s) : \"",paste(missing.var,collapse ="\" \""),"\"\n")
            }

    }
    
    ## ** check covariance links
    if(any(directive==FALSE)){
        if(any(restricted2[directive==FALSE,1] %in% exogenous(object)) || any(restricted2[directive==FALSE,2] %in% exogenous(object))){
            wrong <- union(which(restricted2[directive==FALSE,1] %in% exogenous(object)),
                           which(restricted2[directive==FALSE,2] %in% exogenous(object)))
            stop("covariance link with a exogenous variable: \n",
                 paste(link[wrong], collapse = " ; "),"\n")
        }
    }
    
     return(list(object = object,
                link = link,
                directive = directive,
                restricted = restricted))
}
## * .oneStep_scoresearch
.oneStep_scoresearch  <- function(object, data,
                                  restricted, link, directive,
                                  method.p.adjust, alpha,
                                  cl, trace){

    ## ** initialization
    n.link <- NROW(restricted)
    coef.object <- coef(object)
    namecoef.object <- names(coef.object)
    ncoef.object <- length(coef.object)
    
    ## ** warper
    warper <- function(iterI){ # iterI <- 1

        out <- list(table = data.frame(statistic = as.numeric(NA),
                                       se = as.numeric(NA),
                                       p.value = as.numeric(NA),
                                       adjusted.p.value = as.numeric(NA),
                                       stringsAsFactors = FALSE),
                    iid = NULL)

        ## *** define extended model
        newModel <- addLink(object$model, var1 = restricted[iterI,1], var2 = restricted[iterI,2],
                            covariance = 1-directive[iterI])

        ## *** define constrained coefficients
        coef0.new <- setNames(rep(0, ncoef.object+1), coef(newModel))
        coef0.new[namecoef.object] <- coef.object

        ## *** compute the iid decomposition
        iid.score <- lava::score(newModel, p = coef0.new, data = data, indiv = TRUE)
        Info <- lava::information(newModel, p = coef0.new, n = NROW(data), data = data)

        ## inverse of the information matrix
        sqrt.InfoM1 <- matrixPower(Info, power  = -1/2, symmetric = TRUE, tol = 1e-15)
        ## sqrt.InfoM1 <- solve(chol(Info))

        ## iid decomposition of the normalized score (follows a standard normal distribution)
        iid.normScore <- (iid.score %*% sqrt.InfoM1)
        ## normalized score
        normScore <- colSums(iid.normScore)        

        ## *** store statistics, se, and iid decomposition, and p.value
        out$iid <- rowSums(sweep(iid.normScore, MARGIN = 2, STATS = normScore, FUN = "*"))
        out$iid <- out$iid/sqrt(sum(out$iid^2))
        out$table$statistic <- sqrt(crossprod(normScore))
        out$table$se <- sqrt(sum(out$iid^2))
        return(out)
    }

    ## ** compute score tests
    if(trace>0){
        cat(" - compute score test for all possible additional links \n")
    }
    
    if(!is.null(cl)){
    
        if(trace > 0){
            pb <- utils::txtProgressBar(max = n.link, style = 3) 
        }

        
        ## get influence function
        i <- NULL # [:for CRAN check] foreach
        res <- foreach::`%dopar%`(
                            foreach::foreach(i = 1:n.link,
                                             .combine = function(res1,res2){
                                                 res <- list(df = rbind(res1$table,res2$table),
                                                             iid = cbind(res1$iid,res2$iid))
                                                 return(res)
                                             }), {
                                                 if(trace>0){utils::setTxtProgressBar(pb, i)}
                                                 return(warper(i))
                                             })

        if(trace>0){close(pb)}
        
    }else{
        
        if(trace>0){
            resApply <- pbapply::pblapply(1:n.link, warper)            
        }else{
            resApply <- lapply(1:n.link, warper)
        }
        res <- list(table = do.call(rbind, lapply(resApply,"[[","table")),
                    iid = lapply(resApply,"[[","iid"))
        
    }

    ## index.iid <- unlist(lapply(1:n.link, function(iL){ ## iL <- 1
        ## rep(iL,times = NCOL(res$iid[[iL]]))
    ## }))
    table.test <- data.frame(link = link, res$table, stringsAsFactors = FALSE)
    iid.link <- do.call(cbind,res$iid)

    ## ** p.value
    statistic <- as.numeric(table.test[,"statistic"])
    if(any(statistic<0)){
        stop("Negative score statistic \n")
    }
    ## univariate rejection area
    table.test[,"p.value"] <- 2*(1-stats::pnorm(statistic))

    ## ** adjusted p.value
    if(method.p.adjust == "fastmax"){
        index.max <- which.max(statistic)
        resQmax <- calcDistMaxIntegral(statistic = statistic[index.max[1]], iid = iid.link, df = NULL, alpha = alpha, cl  = cl, trace = trace)

        table.test[index.max, "adjusted.p.value"] <- resQmax$p.adjust
        table.test[index.max, "quantile"] <- resQmax$z
        Sigma <- resQmax$Sigma

    }else if(method.p.adjust == "max"){
        resQmax <- calcDistMaxIntegral(statistic = statistic, iid = iid.link, df = NULL, alpha = alpha, cl  = cl, trace = trace)
        ## resQmax <- calcDistMaxResampling(statistic = statistic, iid = iid.link, index.iid = index.iid,
        ## method.resampling = method.resampling, type.statistic = "chisq",
        ## n.sim = n.sim, alpha = alpha, cl = cl, trace = trace)

        table.test[, "adjusted.p.value"] <- resQmax$p.adjust
        table.test[, "quantile"] <- resQmax$z
        Sigma <- resQmax$Sigma
                
    }else{
        table.test[, "adjusted.p.value"] <- stats::p.adjust(table.test$p.value, method = method.p.adjust)
        table.test[, "quantile"] <- as.numeric(NA)
        Sigma <- NULL        
    }    

    return(list(test = table.test,
                Sigma = Sigma,
                iid = iid.link))
}





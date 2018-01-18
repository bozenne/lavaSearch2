## * Documentation - modelsearch2
#' @title Data-driven Extension of a Latent Variable Model
#' @description Procedure adding relationship between variables that are supported by the data.
#' @name modelsearch2
#' 
#' @param object a \code{lvmfit} object.
#' @param link the name of the additional relationships to consider when expanding the model. Should be a vector containing strings like "Y~X". Optional for \code{lvmfit} objects, see details.
#' @param data [optional] the dataset used to identify the model
#' @param statistic statistic used to perform the test. Can the likelihood ratio test (\code{"LR"}), the score (\code{"score"}) or the max statistic (\code{"max"}).
#' @param method.p.adjust the method used to adjust the p.values for multiple comparisons. Ignored when using the max statistic. Can be any method that is valid for the \code{stats::p.adjust} function (e.g. \code{"fdr"}).
#' @param typeSD [relevant when statistic is Wald] the type of standard error to be used to compute the Wald statistic.
#' Can be \code{"information"}, \code{"robust"} or \code{"jackknife"}.
#' @param df [relevant when statistic is Wald] small sample correction: should the degree of freedom be computed using the Satterthwaite approximation.
#' @param adjust.residuals [relevant when statistic is Wald] small sample correction: should the leverage-adjusted residuals be used to compute the influence function? Otherwise the raw residuals will be used.
#' @param trace should the execution be traced?
#' @param ... additional arguments to be passed to \code{\link{findNewLink}} and \code{.modelsearch2}, see details.
#'
#' @details
#' Argument \code{link}:
#' \itemize{
#' \item \code{lvmfit} object: when not specified all possible additional links are considered.
#' \item other objects: this argument must be specified.
#' }
#'
#' Argument \code{...} passed to \code{\link{findNewLink}}, see the documentation of this function:
#' \itemize{
#' \item exclude.var
#' \item rm.latent_latent
#' \item rm.endo_endo
#' \item rm.latent_endo
#' }
#' 
#' Argument \code{...} passed to \code{\link{modelsearch2}}:
#' \itemize{
#' \item alpha: the significance threshold for retaining a new link.
#' \item method.max: the method used to compute the distribution of the max statistic. See lava.options()$search.calcMaxDist.
#' \item ncpus: the number of cpus that can be used for the computations.
#' \item nStep: the maximum number of links that can be added to the model.
#' \item na.omit: should model leading to NA for the test statistic be ignored. Otherwise this will stop the selection process.
#' }
#' 
#' 
#' @return a latent variable model
#' 
#' @examples
#'
#' #### linear regression ####
#' set.seed(10)
#' mSim <- lvm(Y~X1+X2+X3+X4)
#' addvar(mSim) <- ~Z1+Z2
#' d <- lava::sim(mSim,1e2)
#' eLM <- lm(Y~X1,data = d)
#'
#' res <- modelsearch2(eLM, link = c("Y~X2","Y~X3","Y~X4","Y~Z1","Y~Z2"),
#'              statistic = "LR", method.p.adjust = "holm")
#' res <- modelsearch2(eLM, link = c("Y~X2","Y~X3","Y~X4","Y~Z1","Y~Z2"),
#'              statistic = "Wald", method.p.adjust = "holm", nStep = 1)
#' \dontrun{
#' res <- modelsearch2(eLM, link = c("Y~X2","Y~X3","Y~X4","Y~Z1","Y~Z2"))
#' res <- modelsearch2(eLM, link = c("Y~X2","Y~X3","Y~X4","Y~Z1","Y~Z2"),
#'                     nStep = 2, conditional = TRUE)
#' res$sequenceTest[[2]]
#' }
#' 
#' #### Cox model ####
#' library(survival)
#' data(Melanoma, package = "riskRegression")
#' m <- coxph(Surv(time,status==1)~ici+age, data = Melanoma, x = TRUE, y = TRUE)
#' \dontrun{
#' res <- modelsearch2(m, link = c(status~epicel,status~sex),
#'                     packages = "survival", nStep = 1)
#' }
#' res
#' 
#' #### LVM ####
#' mSim <- lvm()
#' regression(mSim) <- c(y1,y2,y3)~u
#' regression(mSim) <- u~x1+x2
#' categorical(mSim,labels=c("A","B","C")) <- "x2"
#' latent(mSim) <- ~u
#' covariance(mSim) <- y1~y2
#' df <- lava::sim(mSim, 1e2)
#' df$Id <- 1:NROW(df)
#' 
#' m <- lvm(c(y1,y2,y3)~u)
#' latent(m) <- ~u
#' addvar(m) <- ~x1+x2 
#'
#' e <- estimate(m, df)
#' 
#' e$model$index$P1
#'
#' \dontshow{
#' links <- c(u~x1,u~x2C,y3~x2C)
#' resScore <- modelsearch2(e, statistic = "score", link = links, method.p.adjust = "holm")
#' resLR <- modelsearch2(e, statistic = "LR", link = links, method.p.adjust = "holm", nStep = 1)
#' resMax <- modelsearch2(e, rm.endo_endo = TRUE, statistic = "Wald", link = links, nStep = 1)
#' } 
#' \dontrun{
#' resScore <- modelsearch2(e, statistic = "score", method.p.adjust = "holm")
#' resLR <- modelsearch2(e, statistic = "LR", method.p.adjust = "holm")
#' resMax <- modelsearch2(e, rm.endo_endo = TRUE, statistic = "Wald")
#' }
#' 
#'
#' @export
`modelsearch2` <-
  function(object, ...) UseMethod("modelsearch2")

## * method modelsearch2.lvmfit
#' @rdname modelsearch2
#' @export
modelsearch2.lvmfit <- function(object, link = NULL, data = NULL, 
                                statistic = "Wald",  method.p.adjust = "max",
                                typeSD = "information", df = FALSE, adjust.residuals = FALSE,
                                trace = TRUE,
                                ...){

    ## ** normalise arguments
    typeSD <- match.arg(typeSD, c("information","robust","jackknife"))
    if(df == FALSE && adjust.residuals == TRUE){
        stop("Argument \'df\' must be TRUE when arguemnt \'adjust.residuals\' is TRUE \n")
    }

    dots <- list(...)

    args.findNewLink <- c("exclude.var", "rm.latent_latent", "rm.endo_endo", "rm.latent_endo")
    dots.findNewLink <- dots[names(dots) %in% args.findNewLink]

    dots.modelsearch2 <- dots[names(dots) %in% args.findNewLink == FALSE]
    
    ## ** normalize the links
    if(is.null(link)){
        res.find <- do.call(findNewLink,
                            args = c(list(object$model,
                                          data = stats::model.frame(object),
                                          output = "names"),
                                     dots.findNewLink))
        directive <- res.find$directional
        restricted <- res.find$M.links
        link <- res.find$link
        if(is.null(link)){
            stop("Automatic search has not found any possible additional link \n",
                 "Consider specifying manually the argument \'link\' \n")
        }
    }else{
        resLink <- .initializeLinks(object, data = data, link = link,
                                    statistic = statistic)
        object <- resLink$object
        link <- resLink$link
        directive <- resLink$directive
        restricted <- resLink$restricted
    }    
        
    ## ** arguments of the call
    add.args <- setdiff(names(object$call), c("","object","data","control"))
    ls.call <- lapply(add.args, function(arg){object$call[[arg]]})
    names(ls.call) <- add.args

    ls.call$data <- as.data.table(stats::model.frame(object))
    if(!is.null(data)){
        index.cols <- which(names(data)%in%names(ls.call$data)==FALSE)
        if(length(index.cols)>0){
            ls.call$data <- cbind(ls.call$data,
                                  as.data.table(data)[,.SDcols=index.cols])
        }
    }
    ls.call$control <- object$control
    ls.call$control$trace <- FALSE
    
    ## ** extract influence function
    if(typeSD == "jackknife"){
        iid.FCT <- function(x){
            iidJack(x, keep.warnings = FALSE, keep.error = FALSE, trace = FALSE)
        }
        attr(iid.FCT,"method.iid") <- "iidJack"
    }else if(adjust.residuals == FALSE){
        iid.FCT <- function(x){
            res <- lava::iid(x)
            attr(res, "bread") <- NULL
            return(res)
        }
        attr(iid.FCT,"method.iid") <- "iid"
    }else if(typeSD == "iid2"){
        iid.FCT <- function(x){
            iid2(x, adjust.residuals = TRUE)
        }
        attr(iid.FCT,"method.iid") <- "iid2"
    }
    attr(iid.FCT,"typeSD") <- typeSD
    attr(iid.FCT,"df") <- df
    attr(iid.FCT,"adjust.residuals") <- adjust.residuals

    ## ** run modelsearch
    out <- do.call(.modelsearch2,
                   args = c(list(object, 
                                 link = link, restricted = restricted, directive = directive,
                                 statistic = statistic, method.p.adjust = method.p.adjust,
                                 trace = trace,
                                 update.FCT = .updateModelLink.lvm, update.args = ls.call, iid.FCT = iid.FCT),
                            dots.modelsearch2))

    ## ** export
    return(out)
}

## * method modelsearch2.default
#' @rdname modelsearch2
#' @export
modelsearch2.default <- function(object, link, data = NULL,
                                 statistic = "Wald", method.p.adjust = "max", 
                                 typeSD = "information", df = FALSE, adjust.residuals = FALSE,
                                 trace = TRUE,
                                 ...){

    ## ** normalise arguments
    typeSD <- match.arg(typeSD, c("information","robust","jackknife"))
    if(df == FALSE && adjust.residuals == TRUE){
        stop("Argument \'df\' must be TRUE when arguemnt \'adjust.residuals\' is TRUE \n")
    }

    if("lvm" %in% class(object)){
        stop("Cannot apply modelsearch2 to a \"lvm\" object \n",
             "Run the function \'estimate\' first and then \'modelsearch2\' \n")
    }
    if(statistic == "score"){
        stop("Can only use the score statistic with lvmfit objects \n",
             "Consider changing argument \'statistic\' \n")
    }
    if(class(object) %in% c("coxph","cph","phreg")){            
        if(df == TRUE){
            stop("argument \'df\' must be FALSE for Cox models \n")
        }
        if(adjust.residuals == TRUE){
            stop("argument \'adjust.residuals\' must be FALSE for Cox models \n")
        }        
    }else if (!any(paste("score", class(object), sep = ".") %in% methods("score"))) {        
        stop("Extraction of the iid decomposition failed \n",
             "No iid method for models of class ",class(object)," \n")
    }
      
    ## ** get data
    if(is.null(data)){
        data <- evalInParentEnv(object$call$data, envir = environment())
        if(is.null(data)){
            stop("object$call$data not found in the current environment or its parents \n",
                 "consider specify the argument \'data\' \n")
        }
    }

    ## ** get model vars
    model.var <- all.vars(stats::formula(object))
    
    ## ** normalize the links
    restricted <- do.call(cbind,initVarLinks(link, Slink = "~"))
    allVars <- union(model.var, names(data))
    link <- restricted[,2] # does not handle categorical variables    
    directive <- rep(TRUE, length(link))

    if(any(unique(as.vector(restricted)) %in% allVars == FALSE)){
        wrong.var <- unique(as.vector(restricted))[unique(as.vector(restricted)) %in% allVars == FALSE]
        stop("Some links contains variables that are not in the model \n",
             "variables(s) : \"",paste(wrong.var,collapse ="\" \""),"\"\n")
    }
    
    if(any(restricted[,2] %in% model.var)){
        stop("link(s) already in the models: ",paste0(link[restricted[,2] %in% model.var],collapse = " "),"\n")
    }

    ## take care of categorical variables
    for(iLink in 1:length(link)){ # iLink <- 1
        if(is.numeric(data[[link[iLink]]])==FALSE){
            iUlevel <- levels(as.factor(data[[link[iLink]]]))
            if(length(iUlevel)>2){
                stop("Cannot deal with categorical variables with more than 2 levels \n")
            }else{
                link[iLink] <- paste0(link[iLink], iUlevel[2])
            }
        }        
    }
    
    ## ** arguments of the call
    add.args <- setdiff(names(object$call), c("","formula","data"))
    if(length(add.args)>1){
        ls.call <- lapply(add.args, function(arg){object$call[[arg]]})
        names(ls.call) <- add.args        
    }else{
        ls.call <- list()
    }
    ls.call$data <- data

    ## ** extract influence function
    if(class(object) %in% c("coxph","cph")){
        
        if(typeSD == "jackknife"){
            iid.FCT <- function(x){
                iidJack(x, keep.warnings = FALSE, keep.error = FALSE, trace = FALSE)
            }
            attr(iid.FCT, "method.iid") <- "iidJack"        
        }else{
            tryPkg <- requireNamespace("riskRegression")
            if("try-error" %in% class(tryPkg)){
                stop(tryPkg)
            }else if(utils::packageVersion("riskRegression")<="1.4.3"){
                stop("riskRegression version must be > 1.4.3 \n",
                     "latest version available on Github at tagteam/riskRegression \n")
            }else{
                iid.FCT <- function(x){
                    riskRegression::iidCox(x, tau.hazard = 0)$IFbeta
                }
                attr(iid.FCT, "method.iid") <- "iid"
            }
        }
    }else if(adjust.residuals == FALSE){
        iid.FCT <- function(x){
            res <- lava::iid(x)
            attr(res, "bread") <- NULL
            return(res)
        }
        attr(iid.FCT, "method.iid") <- "iid"
    }else{
        iid.FCT <- function(x){
            iid2(x, adjust.residuals = TRUE)
        }
        attr(iid.FCT,"method.iid") <- "iid2"        
    }

    attr(iid.FCT,"typeSD") <- typeSD
    attr(iid.FCT,"df") <- df
    attr(iid.FCT,"adjust.residuals") <- adjust.residuals
    
    ## ** run modelsearch
    out <- .modelsearch2(object, link = link, restricted = restricted, directive = directive,
                         statistic = statistic, method.p.adjust = method.p.adjust,
                         trace = trace,
                         update.FCT = .updateModelLink.default, update.args = ls.call, iid.FCT = iid.FCT, ...)


    ## ** export
    return(out)
}    

## * function .modelsearch2
.modelsearch2 <- function(object, link, restricted, directive,
                          statistic, method.p.adjust,
                          update.FCT, update.args, iid.FCT,
                          nStep = NULL, na.omit = TRUE,
                          alpha = 0.05, method.max = "integration", 
                          ncpus = 1, trace = 1,
                          packages = NULL, conditional = FALSE, exposure = NULL, ## not documented
                          display.warnings = TRUE, export.iid = FALSE ## not documented
                          ){

    adjusted.p.value <- selected <- NULL ## [:for CRAN check] data.table

    ## ** preliminary tests
    method.max <- match.arg(method.max, lava.options()$search.calcMaxDist)
    method.p.adjust <- match.arg(method.p.adjust, lava.options()$search.p.adjust)
    statistic <-  match.arg(statistic, choices = lava.options()$search.statistic)
    
    if(any(exposure %in% names(stats::coef(object)) == FALSE)){
        stop("exposure does not correspond to a coefficient in \'object\' \n")
    }    
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    if(method.p.adjust %in% c("fastmax","max") && statistic != "Wald"){
        stop("Adjustment for multiple testing using the distribution of the max statistic \n",
             "is only available for when specifying statistic=\"Wald\" \n",
             "proposed statistic: ",statistic,"\n",
             "consider changing the value of the argument \'method.p.adjust\' \n")
    }
    if(statistic != "Wald"){
        typeSD <- NULL
    }

    ## ** initialisation
    if(is.null(nStep)){
        nStep <- NROW(restricted)
    }
    iStep <- 1
    iRestricted <- restricted
    iDirective <- directive
    iLink <- link
    iObject <- object
    nObs <- object$data$n
    vec.seqQuantile <- NULL # for conditional testing
    
    ## output
    ls.seqTests <- list()
    ls.seqModels <- list()
    ls.seqIID <- list() # only for method.p.adjust = "max"
    ls.seqSigma <- list() # only for method.p.adjust = "max"

    ## criterion    
    cv <- FALSE
       
    ## cpus
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    if(ncpus>1){
        cl <- parallel::makeCluster(ncpus)
        doParallel::registerDoParallel(cl)
    }
    
    ## ** display a summary of the call
    if(trace>0){
        cat("\n",
            "** Sequential variable selection using the ",statistic," statistic ** \n",
            " Number of possible additional links         : ",length(link)," \n",
            " Maximum number of steps                     : ",nStep,"\n",        
            " Adjustment method for multiple comparisons  : ",method.p.adjust,"\n",
            if(method.p.adjust %in% c("fastmax","max")){
                paste0(" Type of standard error                      : ",attr(iid.FCT,"typeSD"),"\n",
                       " Method for extracting the influence function: ",attr(iid.FCT,"method.iid"),"\n",
                       " Method for computing the quantile           : ",method.max,"\n",
                       " Correction for sequential testing           : ",conditional,"\n")
            },
            " Confidence level                            : ",1-alpha,"\n",
            " Number of cpus                              : ",ncpus,"\n\n",
            sep="")
    }

    ## ** loop over the models
    while(iStep <= nStep && NROW(iRestricted)>0 && cv==FALSE){

        if(trace >= 1){cat("Step ",iStep,":\n",sep="")}        

        if(statistic == "score"){
            ### *** run modelsearch
            res.search <- modelsearch(iObject,
                                      link = gsub("~~","~",iLink),
                                      silent = (trace <= 1))
            
            iN.link <- NROW(iRestricted)        
            res.search$dt.test <- data.table("link" = iLink,
                                             "statistic" = as.numeric(rep(NA,iN.link)),
                                             "p.value" = as.numeric(rep(NA,iN.link)),
                                             "adjusted.p.value" = as.numeric(rep(NA,iN.link)),
                                             "convergence" = 1,
                                             "coefBeta" = as.numeric(rep(NA,iN.link)),
                                             "quantile" = as.numeric(rep(NA,iN.link)),
                                             "corrected.level" = as.numeric(rep(NA,iN.link))
                                             )
            index.match <- match(gsub("~~","~",iLink), res.search$res[,"Index"])
            ## res.search$res[index.match,"Index"]
            res.search$dt.test[, "statistic" := res.search$test[index.match,"Test Statistic"]]
            res.search$dt.test[, "p.value" := res.search$test[index.match,"P-value"]]
            res.search$dt.test[, "adjusted.p.value" := stats::p.adjust(.SD$p.value, method = method.p.adjust)] ## keep .SD for clarity
        }else if(statistic == "LR"){
            ### *** run modelsearchLR
            res.search <- modelsearchLR(iObject, restricted = iRestricted, link = iLink, directive = iDirective,
                                        update.FCT = update.FCT, update.args = update.args,
                                        method.p.adjust = method.p.adjust, display.warnings = display.warnings, trace = trace-1)

        }else if(statistic == "Wald"){
            ### *** run modelsearchMax
            if(conditional && iStep>1){
                iid.previous <- ls.seqIID[[iStep-1]]
                quantile.previous <- vec.seqQuantile[iStep-1]
            }else{
                iid.previous <- NULL
                quantile.previous <- NULL
            }
            res.search <- modelsearchMax(iObject, restricted = iRestricted, link = iLink, directive = iDirective, packages = packages, alpha = alpha,
                                         update.FCT = update.FCT, update.args = update.args, iid.FCT = iid.FCT,
                                         method.p.adjust = method.p.adjust, method.max = method.max,
                                         iid.previous = iid.previous, quantile.previous = quantile.previous,
                                         export.iid = max(conditional,export.iid), trace = trace-1, ncpus = ncpus, initCpus = FALSE)
        }

        ## ** update according the most significant p.value
        ### *** check convergence
        if(na.omit){
            cv <- all(stats::na.omit(res.search$dt.test[["adjusted.p.value"]]) > alpha)
            test.na <- FALSE
        }else{
            cv <- all(res.search$dt.test[["adjusted.p.value"]] > alpha)
            if(is.na(cv)){                    
                cv <- TRUE
                test.na <- TRUE
            }else{
                test.na <- FALSE
            }
        }

        ### *** identify most promising test
        index.rm <- which.max(abs(res.search$dt.test[["statistic"]]))

        ### *** update the output
        res.search$dt.test[,"selected" := .I==index.rm*(1-cv)]
        res.search$dt.test[,"nTests" := .N]
        setkey(res.search$dt.test,statistic)
        rowSelected <- res.search$dt.test[, .I[selected==TRUE]]
		
        ls.seqTests[[iStep]] <- copy(res.search$dt.test)
        if(method.p.adjust %in% "max"){
            ls.seqIID[[iStep]] <- res.search$iid
            ls.seqSigma[[iStep]] <- res.search$Sigma
            vec.seqQuantile <- c(vec.seqQuantile,unique(res.search$dt.test$quantile))
        }
        
        ### *** update the model
        if(cv==FALSE){
            iObject <- update.FCT(iObject, args = update.args,
                                  restricted = iRestricted[index.rm,], directive = iDirective[index.rm])
            if(!is.null(iLink)){iLink <- iLink[-index.rm]}
            iRestricted <- iRestricted[-index.rm,,drop=FALSE]
            iDirective <- iDirective[-index.rm]            
        }
        ls.seqModels[[iStep]] <- iObject

        ### *** display results
        if(trace > 0){
            if(cv==FALSE){
                
                cat("add ",ls.seqTests[[iStep]][rowSelected, link],
                    " (statistic = ",ls.seqTests[[iStep]][rowSelected, statistic],
                    ", adjusted.p.value = ",ls.seqTests[[iStep]][rowSelected, adjusted.p.value],
                    ")\n",sep="")
            }else{
                if(test.na){
                    cat("NA among the test statistics \n")
                }else{
                    rowSelected <- which.max(ls.seqTests[[iStep]][, statistic])
                    cat("no variable to add",
                        " (statistic = ",ls.seqTests[[iStep]][rowSelected, statistic],
                        ", adjusted.p.value = ",ls.seqTests[[iStep]][rowSelected, adjusted.p.value],
                        ")\n",sep="")
                }
            }
        }
        iStep <- iStep + 1
        
    }
         
    ## * Test treatment effect
    if(!is.null(exposure) && FALSE){
        if(statistic == "Wald" && method.p.adjust == "max"){
            ## dt.exposure <- data.table("link" = exposure,
        ##                               "statistic" = numeric(1),
        ##                               "p.value" = numeric(1),
        ##                               "statistic Wald" = numeric(1),
        ##                               "p.value Wald" = numeric(1)
        ##                               )
               
        ##     #### compute statistic
        ##     newBeta <- coef(iObject)[exposure]
        ##     IF.beta <- sqrt(n)*iid(iObject)[,exposure,drop=FALSE]
        ##     SeBeta <- sd(IF.beta)
        ##     dt.exposure[,`statistic` := abs(newBeta/SeBeta)]
            
            
        ##     dt.exposure[,"statistic Wald" := summary(iObject)$coef[exposure,1]/summary(iObject)$coef[exposure,2]]
        ##     dt.exposure[,"p.value Wald" := summary(iObject)$coef[exposure,4]]
            
        ##     #### compute quantile
        ##     if(method.max == "normal"){
        ##         sigma <- cov(iid.conditional)
        ##     }else{
        ##         sigma <- NULL
        ##     }
            
        ##     resQmax <- calcDistMax(dt.exposure[["statistic"]], 
        ##                            iid = cbind(iid.conditional,IF.beta), 
        ##                            mu = c(mu.conditional,newBeta), 
        ##                            conditional = c(mu.conditional,rep(0,length(exposure))),
        ##                            method = method.max,
        ##                            n.sim = n.sim, ncpus = ncpus, initCpus = FALSE, trace = trace)
            
        ##     dt.exposure[,`p.value` := resQmax$p.adjust]
        ##     z <- resQmax$z
            
            ## }else{
            ##   if(identical(iObject$call$robust,TRUE)){
            ##     dt.exposure <- data.table(link = exposure, summary(iObject)$coef[exposure,c("Estimate","Robust SE","P-value"),drop=FALSE])
            ##     dt.exposure[, c("Robust SE") := .SD[["Estimate"]]/.SD[["Robust SE"]], .SDcols = c("Estimate","Robust SE")]
            ##   }else{
            ##     dt.exposure <- data.table(link = exposure, summary(iObject)$coef[exposure,c("Estimate","Z-value","P-value"),drop=FALSE])
            ##   }
            ##     names(dt.exposure) <- c("link","estimate","statistic","p.value")
            ##     z <- NULL
            ## }
        }
    }
    ## * end job
    if(ncpus>1){
        parallel::stopCluster(cl)
    }
    
    ## * export
    if(length(ls.seqIID)==0){ls.seqIID <- NULL}
    if(length(ls.seqSigma)==0){ls.seqSigma <- NULL}
    output <- list(sequenceTest = ls.seqTests,
                   sequenceModel = ls.seqModels,
                   sequenceQuantile = vec.seqQuantile,
                   sequenceIID = ls.seqIID,
                   sequenceSigma = ls.seqSigma,
                   statistic = statistic,
                   method.p.adjust = method.p.adjust,
                   typeSD = attr(iid.FCT,"typeSD"),
                   alpha = alpha,
                   cv = cv)
    class(output) <- "modelsearch2"
    return(output)
}




## * .updateModelLink.lvm
.updateModelLink.lvm <- function(object, args, restricted, directive){
    args$x <- addLink(object$model, var1 = restricted[1], var2 = restricted[2],
                      covariance = 1-directive)

    ## first attempt
    args$start <- stats::coef(object)
    suppressWarnings(
        newObject <- tryCatch(do.call(estimate, args = args),
                         error = function(x){NA},
                         finally = function(x){x})
    )
        
    ## second attempt
    if(newObject$opt$convergence>0){
        args$control$start <- NULL
        suppressWarnings(
            newObject <- tryCatch(do.call(estimate, args = args),
                                  error = function(x){NA},
                                  finally = function(x){x})
        )
    }
    return(newObject)
}

## * .updateModelLink.default
.updateModelLink.default <- function(object, args, restricted, directive, ...){
    FCT.estimate <- as.character(object$call[[1]])

    ## update the formula
    f <- stats::formula(object) #evalInParentEnv(object$call$formula, envir = environment())
    if(is.list(f)){
        test.Y <-  lapply(f, function(ff){
            restricted[1] %in% selectResponse(ff)
        })
        index.Y <- which(unlist(test.Y))
        f[[index.Y]] <- stats::update(f[[index.Y]], stats::as.formula(paste0(".~",restricted[2])))
            
    }else{
        f <- stats::update(f, stats::as.formula(paste0(".~.+",restricted[2])))
    }

    ## update the arguments
    if(length(args)>0){
        for(callArg in names(args)){
            object$call[[callArg]] <- args[[callArg]]
        }
    }
        
    suppressWarnings(
        newObject <- stats::update(object, formula = f)
    )
    newObject$opt$convergence <- 0
    return(newObject)
}



## * .initializeLinks
.initializeLinks <- function(object, data, link, statistic){
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
                   X = var2dummy(object$model, data = data, var = restricted[row,2]), dir = directive[row])
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

            ## modelsearch require to have the variable corresponding to the links in the model
            if(statistic == "score" && any(allVars.data %in% missing.var)){
                addVars <- allVars.data[allVars.data %in% missing.var]
                ff <- stats::as.formula(paste0("~",paste(addVars, collapse = " + ")))
                addvar(object$model) <- ff
                allVars.model <- vars(object$model)        
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
## * Documentation - modelsearch2
#' @title Automatic extension of the lvm
#' @description Add all relevant path between the variables present in a lvm.
#' 
#' @name modelsearch2
#' 
#' @param x a lvm model
#' @param data the dataset used to identify the model
#' @param statistic test statistic used to perform the test. Can the likelihood ratio test (LR), the score (score) or the max statistic (max).
#' @param exposure the links to be tested at the end of the variable selection process. Not functional for now.
#' @param link the links to consider when expanding the model
#' @param exclude.var see the documentation of \code{\link{findNewLink}}.
#' @param rm.latent_latent see the documentation of \code{\link{findNewLink}}.
#' @param rm.endo_endo see the documentation of \code{\link{findNewLink}}.
#' @param rm.latent_endo see the documentation of \code{\link{findNewLink}}.
#' @param nStep the maximum number of links that can be added to the model.
#' @param na.omit should model leading to NA for the test statistic be ignored. Otherwise this will stop the selection process.
#' @param alpha the significance threshold for retaining a new link.
#' @param method.p.adjust the method used to adjust the p.values for multiple comparisons. Ignore when using the max statistic.
#' @param method.max the method used to compute the distribution of the max statistic \code{integration}.
#' @param method.iid the method used to compute the influence function. Can be \code{"iid"} (first order approximation) or \code{"iidJack"} (jacknife estimate).
#' @param conditional should the p.values be corrected for sequential testing. Experimental.
#' @param ncpus the number of cpus that can be used for the computations.
#' @param display.warnings should warnings be display? May occur when dealing with categorical variables or when fitting an extended model.
#' @param trace should the execution be traced?
#' @param packages name of the R packages to be exported to each CPU.
#' @param export.iid should the iid decomposition of the retained coefficient be export. Only relevant when statistic is set to max.
#' @param ... additional arguments to be passed to lower levels functions.
#'
#' @param restricted internal argument.
#' @param directive internal argument.
#' @param update.FCT internal argument.
#' @param update.args internal argument.
#' @param iid.FCT internal argument.
#' 
#' @return a latent variable model
#' 
#' @examples
#'
#' # to save computation time 
#' method.iid <- "iid"
#' 
#' #### linear regression ####
#' mSim <- lvm(Y~X1+X2+X3+X4)
#' addvar(mSim) <- ~Z1+Z2
#' d <- lava::sim(mSim,1e2)
#' eLM <- lm(Y~X1,data = d)
#'
#' res <- modelsearch2(eLM, link = c("Y~X2","Y~X3","Y~X4","Y~Z1","Y~Z2"),
#'              statistic = "LR", method.p.adjust = "holm")
#' res <- modelsearch2(eLM, link = c("Y~X2","Y~X3","Y~X4","Y~Z1","Y~Z2"),
#'              statistic = "Wald", method.p.adjust = "holm",
#'              method.iid = method.iid)
#' \dontrun{
#' res <- modelsearch2(eLM, link = c("Y~X2","Y~X3","Y~X4","Y~Z1","Y~Z2"),
#'                     method.iid = method.iid)
#' res <- modelsearch2(eLM, link = c("Y~X2","Y~X3","Y~X4","Y~Z1","Y~Z2"),
#'                     nStep = 2, 
#'                     method.iid = method.iid, conditional = TRUE)
#' res$sequenceTest[[2]]
#' }
#' 
#' #### Cox model ####
#' library(survival)
#' data(Melanoma, package = "riskRegression")
#' m <- coxph(Surv(time,status==1)~ici+age, data = Melanoma, x = TRUE, y = TRUE)
#' res <- modelsearch2(m, link = c(status~epicel,status~sex),
#'              method.iid = method.iid, packages = "survival")
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
#' ## score
#' \dontrun{
#' resScore <- modelsearch2(e, statistic = "score", method.p.adjust = "holm")
#' resLR <- modelsearch2(e, statistic = "LR", method.p.adjust = "holm")
#' resMax <- modelsearch2(e, rm.endo_endo = TRUE,
#'                         method.iid = method.iid, statistic = "Wald")
#' # resMax <- modelsearch2(e, rm.endo_endo = TRUE,
#' #                        method.iid = method.iid, statistic = "Wald", conditional = TRUE)
#' }
#' \dontshow{
#' links <- c(u~x1,u~x2C,y3~x2C)
#' resScore <- modelsearch2(e, statistic = "score", link = links, method.p.adjust = "holm")
#' resLR <- modelsearch2(e, statistic = "LR", link = links, method.p.adjust = "holm")
#' resMax <- modelsearch2(e, rm.endo_endo = TRUE, statistic = "Wald", method.iid = "iid", link = links)
#' }
#'
#' @export
`modelsearch2` <-
  function(x, ...) UseMethod("modelsearch2")

## * method modelsearch2.lvmfit
#' @rdname modelsearch2
#' @export
modelsearch2.lvmfit <- function(x, data = NULL, link = NULL,
                                method.iid = "iid",
                                exclude.var = NULL, rm.latent_latent= FALSE, rm.endo_endo= FALSE, rm.latent_endo= FALSE,
                                ...){

    ## ** normalise arguments
    method.iid <- match.arg(method.iid, lava.options()$search.iid)

    ## ** normalize the links
    if(is.null(link)){
        res.find <- findNewLink(x$model,
                                data = model.frame(x),
                                exclude.var = exclude.var,
                                rm.latent_latent = rm.latent_latent,
                                rm.endo_endo = rm.endo_endo, rm.latent_endo = rm.latent_endo,
                                output = "names")
        directive <- res.find$directional
        restricted <- res.find$M.links
        link <- res.find$link
        if(is.null(link)){
            stop("Automatic search has not found any possible additional link \n",
                 "Consider specifying manually the argument \'link\' \n")
        }
    }else{
        restricted <- do.call(cbind,initVarLinks(link))
        directive <- rep(TRUE, length(link))
        index.Ndir <- grep(lava.options()$symbols[2],link,fixed=TRUE)
        if(length(index.Ndir)>0){
            directive[index.Ndir] <- FALSE
        }

        ## all vars
        if(is.null(data)){            
            data <- evalInParentEnv(x$call$data, envir = environment())             
            if(is.null(data)){
                data <- lava::sim(x,1)
            }            
        }
        ## take care of categorical variables
        ls.linkvar <- do.call(rbind,lapply(1:NROW(restricted), function(row){
            data.frame(Y = restricted[row,1],
                       X = var2dummy(x$model, data = data, var = restricted[row,2]), dir = directive[row])
        }))
        restricted2 <- as.matrix(ls.linkvar[,1:2,drop=FALSE])
        directive <- ls.linkvar[,3]
        link <- paste0(restricted2[,1],lava.options()$symbols[2-directive],restricted2[,2])

        ## check links
        allVars.link <- unique(as.vector(restricted2))
        allVars.model <- vars(x$model)
        allVars.data <- names(data)
        if(any(allVars.link %in% allVars.model == FALSE)){
            missing.var <- allVars.link[allVars.link %in% allVars.model == FALSE]
            
            if(any(allVars.link %in% allVars.data == FALSE)){
                missing.var <- allVars.link[allVars.link %in% allVars.data == FALSE]
                stop("Some links contains variables that are not in the latent variable model \n",
                     "variables(s) : \"",paste(missing.var,collapse ="\" \""),"\"\n")
            }

            ## modelsearch require to have the variable corresponding to the links in the model
            if(list(...)$statistic == "score" && any(allVars.data %in% missing.var)){
                addVars <- allVars.data[allVars.data %in% missing.var]
                ff <- as.formula(paste0("~",paste(addVars, collapse = " + ")))
                addvar(x$model) <- ff
                allVars.model <- vars(x$model)        
            }

            
        }
        
        ## check covariance links
        if(any(directive==FALSE)){
            if(any(restricted2[directive==FALSE,1] %in% exogenous(x)) || any(restricted2[directive==FALSE,2] %in% exogenous(x))){
                wrong <- union(which(restricted2[directive==FALSE,1] %in% exogenous(x)),
                               which(restricted2[directive==FALSE,2] %in% exogenous(x)))
                stop("covariance link with a exogenous variable: \n",
                     paste(link[wrong], collapse = " ; "),"\n")
            }
        }

    }    

    ## ** arguments of the call
    add.args <- setdiff(names(x$call), c("","x","data","control"))
    ls.call <- lapply(add.args, function(arg){x$call[[arg]]})
    names(ls.call) <- add.args

    ls.call$data <- x$data$model.frame
    if(!is.null(data)){
        index.cols <- which(names(data)%in%names(ls.call$data)==FALSE)
        if(length(index.cols)>0){
            ls.call$data <- cbind(ls.call$data,
                                  as.list(data)[index.cols])
        }
    }
    ls.call$control <- x$control
    ls.call$control$trace <- FALSE

    ## ** extract influence function
    if(method.iid == "iid"){
        iid.FCT <- function(x){
            res <- lava::iid(x)
            attr(res, "bread") <- NULL
            return(res)
        }
        attr(iid.FCT,"method.iid") <- "iid"
    }else if(method.iid == "iidJack"){
        iid.FCT <- function(x){
            iidJack(x, keep.warnings = FALSE, keep.error = FALSE, trace = FALSE)
        }
        attr(iid.FCT,"method.iid") <- "iidJack"
    }

    ## ** run modelsearch
    out <- .modelsearch2(x, 
                         link = link, restricted = restricted, directive = directive,
                         update.FCT = .updateModelLink.lvm, update.args = ls.call, iid.FCT = iid.FCT,...)

    ## ** export
    return(out)
}

## * method modelsearch2.default
#' @rdname modelsearch2
#' @export
modelsearch2.default <- function(x, link, data = NULL,
                                 method.iid = "iid", statistic = "Wald", ...){

    ## ** normalise arguments
    method.iid <- match.arg(method.iid, lava.options()$search.iid)
    if("lvm" %in% class(x)){
        stop("Cannot apply modelsearch2 to a \"lvm\" object \n",
             "Run the function \'estimate\' first and then \'modelsearch2\' \n")
    }
    if("lvmfit" %in% class(x) == FALSE && statistic == "score"){
        stop("Can only use the score statistic with lvmfit objects \n",
             "Consider changing argument \'statistic\' \n")
    }
    if(class(x) %in% c("coxph","cph","phreg")){            
        ##
    }else if (!any(paste("score", class(x), sep = ".") %in% methods("score"))) {        
        stop("Extraction of the iid decomposition failed \n",
             "No iid method for models of class ",class(x)," \n")
    }
      
    ## ** get data
    if(is.null(data)){
        data <- evalInParentEnv(x$call$data, envir = environment())
        if(is.null(data)){
            stop("x$call$data not found in the current environment or its parents \n",
                 "consider specify the argument \'data\' \n")
        }
    }

    ## ** get model vars
    model.var <- all.vars(formula(x))
    
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
    add.args <- setdiff(names(x$call), c("","formula","data"))
    if(length(add.args)>1){
        ls.call <- lapply(add.args, function(arg){x$call[[arg]]})
        names(ls.call) <- add.args        
    }else{
        ls.call <- list()
    }
    ls.call$data <- data

    ## ** extract influence function
    if(class(x) %in% c("coxph","cph")){
        
        if(method.iid == "iid"){
            iid.FCT <- function(x){
                riskRegression::iidCox(x, tau.hazard = 0)$IFbeta
            }
            attr(iid.FCT, "method.iid") <- "iid"
        }else if(method.iid == "iidJack"){
            iid.FCT <- function(x){
                iidJack(x, keep.warnings = FALSE, keep.error = FALSE, trace = FALSE)
            }
            attr(iid.FCT, "method.iid") <- "iidJack"
        }
        
    }else{

        if(method.iid == "iid"){
            iid.FCT <- function(x){
                res <- lava::iid(x)
                attr(res, "bread") <- NULL
                return(res)
            }
            attr(iid.FCT, "method.iid") <- "iid"
        }else if(method.iid == "iidJack"){
            iid.FCT <- function(x){
                iidJack(x, keep.warnings = FALSE, keep.error = FALSE, trace = FALSE)
            }
            attr(iid.FCT, "method.iid") <- "iidJack"
        }
        
    } 
    
    ## ** run modelsearch
    out <- .modelsearch2(x, statistic = statistic,
                         link = link, restricted = restricted, directive = directive,
                         update.FCT = .updateModelLink.default, update.args = ls.call, iid.FCT = iid.FCT, ...)


    ## ** export
    return(out)
}    

## * function .modelsearch2
#' @rdname modelsearch2
#' @export
.modelsearch2 <- function(x, link, restricted, directive, update.FCT, update.args, iid.FCT, packages = NULL,
                          statistic = "Wald", conditional = FALSE, exposure = NULL, method.p.adjust = "max",                          
                          nStep = NULL, na.omit = TRUE,
                          alpha = 0.05, method.max = "integration", 
                          ncpus = 1,
                          display.warnings = TRUE, trace = 1, export.iid = FALSE){

    ## ** preliminary tests    
    method.max <- match.arg(method.max, lava.options()$search.calcMaxDist)
    method.p.adjust <- match.arg(method.p.adjust, lava.options()$search.p.adjust)
    statistic <-  match.arg(statistic, choices = lava.options()$search.statistic)

    if(any(exposure %in% names(coef(x)) == FALSE)){
        stop("exposure does not correspond to a coefficient in \'x\' \n")
    }    
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    if(method.p.adjust == "max" && statistic != "Wald"){
        stop("Adjustment for multiple testing using the distribution of the max statistic \n",
             "is only available for when specifying statistic=\"Wald\" \n",
             "proposed statistic: ",statistic,"\n",
             "consider changing the value of the argument \'method.p.adjust\' \n")
    }
    if(statistic != "Wald"){
        method.iid <- NULL
        iid.FCT <- NULL
    }else{
        method.iid <- attr(iid.FCT,"method.iid")
        attr(iid.FCT,"method.iid") <- NULL
    }

    ## ** initialisation
    if(is.null(nStep)){
        nStep <- NROW(restricted)
    }
    iStep <- 1
    iRestricted <- restricted
    iDirective <- directive
    iLink <- link
    iObject <- x
    nObs <- x$data$n
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
        doSNOW::registerDoSNOW(cl)
    }
    
    ## ** display a summary of the call
    if(trace>0){
        cat("\n",
            "** Sequential variable selection using the ",statistic," statistic ** \n",
            " Number of possible additional links         : ",length(link)," \n",
            " Maximum number of steps                     : ",nStep,"\n",        
            " Adjustment method for multiple comparisons  : ",method.p.adjust,"\n",
            if(method.p.adjust=="max"){
                paste0(" Method for extracting the influence function: ",method.iid,"\n",
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
            ## *** run modelsearch
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
            # res.search$res[index.match,"Index"]
            res.search$dt.test[, c("statistic") := res.search$test[index.match,"Test Statistic"]]
            res.search$dt.test[, c("p.value") := res.search$test[index.match,"P-value"]]
            res.search$dt.test[, c("adjusted.p.value") := p.adjust(.SD$p.value, method = method.p.adjust)]
        }else if(statistic == "LR"){
            ## *** run modelsearchLR
            res.search <- modelsearchLR(iObject, restricted = iRestricted, link = iLink, directive = iDirective,
                                        update.FCT = update.FCT, update.args = update.args,
                                        method.p.adjust = method.p.adjust, display.warnings = display.warnings, trace = trace-1)

        }else if(statistic == "Wald"){
            ## *** run modelsearchMax
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
            cv <- all(na.omit(res.search$dt.test[["adjusted.p.value"]]) > alpha)
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

        ### *** identify most pryomising test
        index.rm <- which.max(abs(res.search$dt.test[["statistic"]]))

        ### *** update the output
        res.search$dt.test[,c("selected") := .I==index.rm*(1-cv)]
        res.search$dt.test[,c("nTests") := .N]
        setkey(res.search$dt.test,statistic)
        rowSelected <- res.search$dt.test[, .I[.SD$selected==TRUE]]
		
        ls.seqTests[[iStep]] <- copy(res.search$dt.test)
        if(method.p.adjust == "max"){
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
                
                cat("add ",ls.seqTests[[iStep]][rowSelected, .SD$link],
                    " (statistic = ",ls.seqTests[[iStep]][rowSelected,.SD$statistic],
                    ", adjusted.p.value = ",ls.seqTests[[iStep]][rowSelected,.SD$adjusted.p.value],
                    ")\n",sep="")
            }else{
                if(test.na){
                    cat("NA among the test statistics \n")
                }else{
                    cat("no variable to add \n")
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
                   alpha = alpha,
                   method.iid = method.iid,
                   cv = cv)
    class(output) <- "modelsearch2"
    return(output)
}




## * .updateModelLink.lvm
.updateModelLink.lvm <- function(x, args, restricted, directive){
    args$x <- addLink(x$model, var1 = restricted[1], var2 = restricted[2],
                      covariance = 1-directive)

    ## first attempt
    args$start <- coef(x)
    suppressWarnings(
        newx <- tryCatch(do.call(estimate, args = args),
                         error = function(x){NA},
                         finally = function(x){x})
    )
        
    ## second attempt
    if(newx$opt$convergence>0){
        args$control$start <- NULL
        suppressWarnings(
            newx <- tryCatch(do.call(estimate, args = args),
                             error = function(x){NA},
                             finally = function(x){x})
        )
    }
    return(newx)
}

## * .updateModelLink.default
.updateModelLink.default <- function(x, args, restricted, directive, ...){
    FCT.estimate <- as.character(x$call[[1]])

    ## update the formula
    f <- formula(x) #evalInParentEnv(x$call$formula, envir = environment())
    if(is.list(f)){
        test.Y <-  lapply(f, function(ff){
            restricted[1] %in% selectResponse(ff)
        })
        index.Y <- which(unlist(test.Y))
        f[[index.Y]] <- update(f[[index.Y]], as.formula(paste0(".~",restricted[2])))
            
    }else{
        f <- update(f, as.formula(paste0(".~.+",restricted[2])))
    }

    ## update the arguments
    if(length(args)>0){
        for(callArg in names(args)){
            x$call[[callArg]] <- args[[callArg]]
        }
    }
        
    suppressWarnings(
        newx <- update(x, formula = f)
    )
    newx$opt$convergence <- 0
    return(newx)
}



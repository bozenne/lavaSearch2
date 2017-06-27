# {{{ documentation
#' @title Automatic extension of the lvm
#' @description Add all relevant path between the variables present in a lvm.
#' 
#' @name modelsearch2
#' 
#' @param x a lvm model
#' @param data the dataset used to identify the model
#' @param statistic test statistic used to perform the test. Can the likelihood ratio test (LR), the score (score) or the max statistic (max).
#' @param exposure the links to be tested at the end of the variable selection process.
#' @param link the links to consider when expanding the model
#' @param exclude.var see the documentation of \code{\link{findNewLink}}.
#' @param rm.latent_latent see the documentation of \code{\link{findNewLink}}.
#' @param rm.endo_endo see the documentation of \code{\link{findNewLink}}.
#' @param rm.latent_endo see the documentation of \code{\link{findNewLink}}.
#' @param nStep the maximum number of links that can be added to the model.
#' @param na.omit should model leading to NA for the test statistic be ignored. Otherwise this will stop the selection process.
#' @param alpha the significance threshold for retaining a new link
#' @param method.p.adjust the method used to adjust the p.values for multiple comparisons. Ignore when using the max statistic.
#' @param method.max the method used to compute the distribution of the max statistic \code{integration}.
#' @param method.iid the method used to compute the influence function. Can be \code{"iid"} (first order approximation) or \code{"iidJack"} (jacknife estimate).
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
#' d <- sim(mSim,1e2)
#' eLM <- lm(Y~X1,data = d)
#'
#' res <- modelsearch2(eLM, link = c("Y~X2","Y~X3","Y~X4","Y~Z1","Y~Z2"),
#'              statistic = "LR")
#' res <- modelsearch2(eLM, link = c("Y~X2","Y~X3","Y~X4","Y~Z1","Y~Z2"),
#'              statistic = "max", method.iid = method.iid)
#'
#' #### Cox model ####
#' library(survival)
#' data(Melanoma, package = "riskRegression")
#' m <- coxph(Surv(time,status==1)~ici+age, data = Melanoma, x = TRUE, y = TRUE)
#' modelsearch2(m, link = c(status~epicel,status~sex),
#'              method.iid = method.iid, packages = "survival")
#' 
#' #### LVM ####
#' mSim <- lvm()
#' regression(mSim) <- c(y1,y2,y3)~u
#' regression(mSim) <- u~x1+x2
#' categorical(mSim,labels=c("A","B","C")) <- "x2"
#' latent(mSim) <- ~u
#' covariance(mSim) <- y1~y2
#' df <- sim(mSim, 1e2)
#' df$Id <- 1:NROW(df)
#' 
#' m <- lvm(c(y1,y2,y3)~u)
#' latent(m) <- ~u
#' addvar(m) <- ~x1+x2 
#'
#' e <- estimate(m, df)
#'
#' ## score
#' \dontrun{
#' resScore <- modelsearch2(e, statistic = "score")
#' resLR <- modelsearch2(e, statistic = "LR")
#' resMax <- modelsearch2(e, rm.endo_endo = TRUE, statistic = "max")
#' }
#' \dontshow{
#' links <- c(u~x1,u~x2C,y3~x2C)
#' resScore <- modelsearch2(e, statistic = "score", link = links)
#' resLR <- modelsearch2(e, statistic = "LR", link = links)
#' resMax <- modelsearch2(e, rm.endo_endo = TRUE, statistic = "max", method.iid = "iid", link = links)
#' }
#'
#' @export
`modelsearch2` <-
  function(x, ...) UseMethod("modelsearch2")
# }}}

# {{{ modelsearch2.lvmfit
#' @rdname modelsearch2
#' @export
modelsearch2.lvmfit <- function(x, data = NULL, link = NULL,
                                method.iid = "iidJack",
                                exclude.var = NULL, rm.latent_latent= FALSE, rm.endo_endo= FALSE, rm.latent_endo= FALSE,
                                ...){

    #### normalize the links ####
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
    }else{
        restricted <- do.call(cbind,initVar_links(link))
        
        directive <- rep(TRUE, length(link))
        index.Ndir <- grep(lava.options()$symbols[2],link,fixed=TRUE)
        if(length(index.Ndir)){
            directive[index.Ndir] <- FALSE
        }
        
        ## all vars
        if(is.null(data)){
            allVars <- union(vars(x),vars(lava_categorical2dummy(x, sim(x,1))$x))
        }else{
            allVars <- union(vars(x),vars(lava_categorical2dummy(x, data)$x))
        }
        if(any(unique(as.vector(restricted)) %in% allVars == FALSE)){
            wrong.var <- unique(as.vector(restricted))[unique(as.vector(restricted)) %in% allVars == FALSE]
            stop("Some links contains variables that are not in the latent variable model \n",
                 "variables(s) : \"",paste(wrong.var,collapse ="\" \""),"\"\n")
        }
        ## take care of categorical variables
        ls.linkvar <- do.call(rbind,lapply(1:NROW(restricted), function(row){
            data.frame(Y = restricted[row,1], X = var2dummy(x$model0, data = data, var = restricted[row,2]), dir = directive[row])
        }))
        restricted <- as.matrix(ls.linkvar[,1:2,drop=FALSE])
        directive <- ls.linkvar[,3]
        link <- paste0(restricted[,1],lava.options()$symbols[2-directive],restricted[,2])

        ## check covariance links
        if(any(directive==FALSE)){
            if(any(restricted[directive==FALSE,1] %in% exogenous(x)) || any(restricted[directive==FALSE,2] %in% exogenous(x))){
                wrong <- union(which(restricted[directive==FALSE,1] %in% exogenous(x)),
                               which(restricted[directive==FALSE,2] %in% exogenous(x)))
                stop("covariance link with a exogenous variable: \n",
                     paste(link[wrong], collapse = " ; "),"\n")
            }
        }

    }    

    #### arguments of the call ####
    add.args <- setdiff(names(x$call), c("","x","data","control"))
    ls.call <- lapply(add.args, function(arg){x$call[[arg]]})
    names(ls.call) <- add.args
    if(is.null(data)){
        ls.call$data <- x$data$model.frame
    }else{
        ls.call$data <- data
    }
    ls.call$control <- x$control
    ls.call$control$trace <- FALSE


    #### update the model ####
    update.FCT <- function(x, args, restricted, directive){
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
    # res <- update.x(x, ls.call, restricted = restricted[1,], directive[1])

     iid.FCT <- function(x){
                res <- lava::iid(x)
                attr(res, "bread") <- NULL
                return(res)
     }
    
    #### run modelsearch ####
    out <- .modelsearch2(x, 
                         link = link, restricted = restricted, directive = directive,
                         update.FCT = update.FCT, update.args = ls.call, iid.FCT = iid.FCT,...)

    #### export ####
    return(out)
}
# }}}

# {{{ modelsearch2.default
#' @rdname modelsearch2
#' @export
modelsearch2.default <- function(x, link, data = NULL,
                                 method.iid = "iidJack", statistic = "max", ...){

    if(class(x) %in% c("coxph","cph","phreg")){            
        test.cox <- try(riskRegression::iidCox(x), silent = TRUE)
        if(class(test.cox) == "try-error"){
           print(test.cox)
        }
    }else if (!any(paste("score", class(x), sep = ".") %in% methods("score"))) {        
        stop("Extraction of the iid decomposition failed \n",
                "No iid method for models of class ",class(x)," \n")
    }
    statistic <-  match.arg(statistic, c("max","LR"))
    
    #### get data ####
    if(is.null(data)){
        data <- eval(x$call$data)
    }

    #### get model vars ####
    model.var <- all.vars(formula(x))
    
    #### normalize the links ####
    restricted <- do.call(cbind,initVar_links(link))
    allVars <- union(model.var, names(data))
    link <- restricted[,2] # does not handle categorical variables    
    directive <- TRUE        
    
    if(any(unique(as.vector(restricted)) %in% allVars == FALSE)){
        wrong.var <- unique(as.vector(restricted))[unique(as.vector(restricted)) %in% allVars == FALSE]
        stop("Some links contains variables that are not in the latent variable model \n",
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
    
    #### arguments of the call ####
    add.args <- setdiff(names(x$call), c("","formula","data"))
    if(length(add.args)>1){
        ls.call <- lapply(add.args, function(arg){x$call[[arg]]})
        names(ls.call) <- add.args        
    }else{
        ls.call <- list()
    }
    ls.call$data <- data

    #### update the model ####
    update.FCT <- function(x, args, restricted, directive){
        FCT.estimate <- as.character(x$call[[1]])

        ## update the formula
        f <- eval(x$call$formula)
        if(is.list(f)){
            test.Y <-  lapply(f, function(ff){
                restricted[1] %in% select.response(ff)
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

    #### get iid ####
    if(class(x) %in% c("coxph","cph")){
        
        if(method.iid == "iid"){
            iid.FCT <- function(x){
                riskRegression::iidCox(x, tauHazard = 0)$IFbeta
            }
        }else{
            iid.FCT <- function(x){
                iidJack(x)
            }
        }
        
    }else{

        if(method.iid == "iid"){
            iid.FCT <- function(x){
                res <- lava::iid(x)
                attr(res, "bread") <- NULL
                return(res)
            }
        }else{
            iid.FCT <- function(x){
                iidJack(x)
            }
        }
        
    }
    # iid.FCT(x) 
    
    #### run modelsearch ####
    out <- .modelsearch2(x, statistic = statistic,
                         link = link, restricted = restricted, directive = directive,
                         update.FCT = update.FCT, update.args = ls.call, iid.FCT = iid.FCT, ...)


    #### export ####
    return(out)
}    
# }}}

# {{{ .modelsearch2
#' @rdname modelsearch2
#' @export
.modelsearch2 <- function(x, link, restricted, directive, update.FCT, update.args, iid.FCT, packages = NULL,
                          statistic = "score", exposure = NULL,                           
                          nStep = NULL, na.omit = TRUE,
                          alpha = 0.05, method.max = "integration", method.p.adjust = "bonferroni",
                          ncpus = 1,
                          display.warnings = TRUE, trace = 2, export.iid = FALSE){
    
    # for CRAN check
    `Estimate` <- `Robust SE` <- n.sim <- n <- coefBeta <- p.value <- adjusted.p.value <- NULL
    
    # {{{ preliminary tests    
    method.max <- match.arg(method.max, lava.options()$calcDistMax.method)
    statistic <-  match.arg(statistic, choices = c("score","LR", "max"))
    if(any(exposure %in% names(coef(x)) == FALSE)){
        stop("exposure does not correspond to a coefficient in \'x\' \n")
    }    
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    # }}}

    # {{{ initialisation
    if(is.null(nStep)){
        nStep <- NROW(restricted)
    }
    iStep <- 1
    iRestricted <- restricted
    iDirective <- directive
    iLink <- link
    iObject <- x
    nObs <- x$data$n

    conditional <- NULL # compatibility for all stat when exporting
    iid.conditional <- NULL # compatibility for all stat when exporting
    if(statistic == "max"){        
        mu.conditional <- NULL        
    }

    ## output
    dt.testOverModel <- NULL

    ## criterion    
    cv <- FALSE

    ## cpus
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    cl <- parallel::makeCluster(ncpus)
    doSNOW::registerDoSNOW(cl)
    # }}}
    # {{{ loop over the models
    while(iStep <= nStep && NROW(iRestricted)>0 && cv==FALSE){

        if(trace >= 1){cat("Step ",iStep,":\n",sep="")}        

        if(statistic == "score"){
            # {{{ run modelsearch
            res.search <- modelsearch(iObject, link = iLink, silent = (trace <= 1))

            iN.link <- NROW(iRestricted)        
            res.search$dt.test <- data.table("link" = iLink,
                                             "statistic" = as.numeric(rep(NA,iN.link)),
                                             "p.value" = as.numeric(rep(NA,iN.link)),
                                             "adjusted.p.value" = as.numeric(rep(NA,iN.link))
                                             )
            index.match <- match(gsub("~~","~",iLink), res.search$res[,"Index"])
            # res.search$res[index.match,"Index"]
            res.search$dt.test[, statistic := res.search$test[index.match,"Test Statistic"]]
            res.search$dt.test[, p.value := res.search$test[index.match,"P-value"]]
            res.search$dt.test[, adjusted.p.value := p.adjust(p.value, method = method.p.adjust)]
            # }}}                
        }else if(statistic == "LR"){
            # {{{ run modelsearchLR            
            res.search <- modelsearchLR(iObject, restricted = iRestricted, link = iLink, directive = iDirective,
                                        update.FCT = update.FCT, update.args = update.args,
                                        method.p.adjust = method.p.adjust, display.warnings = display.warnings, trace = trace-1)
    
            # }}}
        }else if(statistic == "max"){
            # {{{ run modelsearchMax
            res.search <- modelsearchMax(iObject, restricted = iRestricted, link = iLink, directive = iDirective, packages = packages, alpha = alpha,
                                         update.FCT = update.FCT, update.args = update.args, iid.FCT = iid.FCT,
                                         method.p.adjust = method.p.adjust, method.max = method.max,
                                         conditional = conditional, mu.conditional = mu.conditional, iid.conditional = iid.conditional,
                                         export.iid = max(1,export.iid), trace = trace-1, ncpus = ncpus, initCpus = FALSE)            # }}}
        }

#        names(res.search)
 #       res.search$dt.test
        # {{{ update according the most significant p.value

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
        
        if(cv==FALSE){            
            index.rm <- which.max(abs(res.search$dt.test[["statistic"]]))

            if(statistic == "max"){
                mu.conditional <- c(mu.conditional, res.search$dt.test[index.rm,coefBeta])
                conditional <- c(conditional,setNames(res.search$z,iLink[index.rm]))        
                iid.conditional <- cbind(iid.conditional, res.search$iid.link)
            }
            ## update the model
            iObject <- update.FCT(iObject, args = update.args,
                                  restricted = iRestricted[index.rm,], directive = iDirective[index.rm])
            if(!is.null(iLink)){iLink <- iLink[-index.rm]}
            iRestricted <- iRestricted[-index.rm,,drop=FALSE]
            iDirective <- iDirective[-index.rm]
            dt.testOverModel <- rbind(dt.testOverModel, res.search$dt.test[index.rm])
        }else{
            if(statistic == "max"){
                conditional <- c(conditional,res.search$z)
            }
        }
        
        if(trace > 0){
            if(cv==FALSE){
                cat("add ",dt.testOverModel[.N,link],"\n")
            }else{
                if(test.na){
                    cat("NA among the test statistics \n")
                }else{
                    cat("no variable to add \n")
                }
            }
        }
        # }}}
        iStep <- iStep + 1

    }
    # }}}
    
    # {{{ Test treatment effect
    if(!is.null(exposure)){
        if(statistic == "max" && length(mu.conditional)>0){
            dt.exposure <- data.table("link" = exposure,
                                      "statistic" = numeric(1),
                                      "p.value" = numeric(1),
                                      "statistic Wald" = numeric(1),
                                      "p.value Wald" = numeric(1)
                                      )
               
            #### compute statistic
            newBeta <- coef(iObject)[exposure]
            IF.beta <- sqrt(n)*iid(iObject)[,exposure,drop=FALSE]
            SeBeta <- sd(IF.beta)
            dt.exposure[,`statistic` := abs(newBeta/SeBeta)]
            
            
            dt.exposure[,"statistic Wald" := summary(iObject)$coef[exposure,1]/summary(iObject)$coef[exposure,2]]
            dt.exposure[,"p.value Wald" := summary(iObject)$coef[exposure,4]]
            
            #### compute quantile
            if(method.max == "normal"){
                sigma <- cov(iid.conditional)
            }else{
                sigma <- NULL
            }
            
            resQmax <- calcDistMax(dt.exposure[["statistic"]], 
                                   iid = cbind(iid.conditional,IF.beta), 
                                   mu = c(mu.conditional,newBeta), 
                                   conditional = c(mu.conditional,rep(0,length(exposure))),
                                   method = method.max,
                                   n.sim = n.sim, ncpus = ncpus, initCpus = FALSE, trace = trace)
            
            dt.exposure[,`p.value` := resQmax$p.adjust]
            z <- resQmax$z
            
        }else{
          if(identical(iObject$call$robust,TRUE)){
            dt.exposure <- data.table(link = exposure, summary(iObject)$coef[exposure,c("Estimate","Robust SE","P-value"),drop=FALSE])
            dt.exposure[, `Robust SE` := `Estimate`/`Robust SE`]
          }else{
            dt.exposure <- data.table(link = exposure, summary(iObject)$coef[exposure,c("Estimate","Z-value","P-value"),drop=FALSE])
          }
            names(dt.exposure) <- c("link","estimate","statistic","p.value")
            z <- NULL
        }
        
        
        
    }else{
        dt.exposure <- NULL
        z <- NULL
    }
    # }}}

    ## end job
    parallel::stopCluster(cl)

    #### export
    output <- list(lvm = iObject,
                   sequenceTest = dt.testOverModel,
                   exposure = dt.exposure,
                   quantile.exposure = z,
                   lastTest = res.search$dt.test,
                   conditional = conditional,
                   iid.conditional = if(export.iid>0){iid.conditional}else{NULL},
                   iid.lastTest = if(export.iid>1){res.search$iid.all}else{NULL},
                   statistic = statistic,
                   cv = cv)
    class(output) <- "modelsearch2"
    return(output)
}


# {{{ print.modelsearch2
#' @method print modelsearch2
#' @export
print.modelsearch2 <- function(x, ...){
    . <- NULL
    
    link <- statistic <- adjusted.p.value <- NULL
    
  if(is.null(x$sequenceTest)){
    cat("The variable selection procedure did not retain any variable \n") 
    print(x$lastTest)
  }else{
    n.var <- NROW(x$sequenceTest)
    cat("The variable selection procedure retained ",n.var," variable",if(n.var>1){"s"},"\n") 
    print(x$sequenceTest[,.(link,statistic,adjusted.p.value)])
  }
    if(x$statistic == "max"){
        cat("significance quantile for the max statistic : ",x$conditional,"\n")
    }
}
# }}}


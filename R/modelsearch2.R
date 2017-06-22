# {{{ modelsearch2
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
#' @param ncpus the number of cpus that can be used for the computations.
#' @param display.warnings should warnings be display? May occur when dealing with categorical variables or when fitting an extended model.
#' @param trace should the execution be traced?
#' @param export.iid should the iid decomposition of the retained coefficient be export. Only relevant when statistic is set to max.
#' @param ... additional arguments to be passed to lower levels functions.
#'
#' @return a latent variable model
#' 
#' @examples 
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
#'
#' e <- estimate(m, df, robust = TRUE, cluster = "Id")
#' resMax <- modelsearch2(e, rm.endo_endo = TRUE, statistic = "max")
#' }
#' \dontshow{
#' links <- c(u~x1,u~x2C,y3~x2C)
#' resScore <- modelsearch2(e, statistic = "score", link = links)
#' resLR <- modelsearch2(e, statistic = "LR", link = links)
#' e <- estimate(m, df, robust = TRUE, cluster = "Id")
#' resMax <- modelsearch2(e, rm.endo_endo = TRUE, statistic = "max", link = links)
#' }
#' set.seed(10)
#' system.time(
#' resMax1 <- modelsearch2(e, rm.endo_endo = TRUE, statistic = "max")
#' )
#' resMax1$lastTest
#' resMax1
#'
#' @export
`modelsearch2` <-
  function(x, ...) UseMethod("modelsearch2")

#' @rdname modelsearch2
#' @export
modelsearch2.lvmfit <- function(x, data = NULL, statistic = "score", exposure = NULL, 
                                link = NULL, exclude.var = NULL, rm.latent_latent= FALSE, rm.endo_endo= FALSE, rm.latent_endo= FALSE,
                                nStep = NULL, na.omit = TRUE,
                                alpha = 0.05, method.max = "integration", method.p.adjust = "bonferroni",
                                ncpus = 1,
                                display.warnings = TRUE, trace = 2, export.iid = FALSE, ...){

    `Estimate` <- `Robust SE` <- n.sim <- n <- coefBeta <- p.value <- adjusted.p.value <- NULL # for CRAN check
    
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    match.arg(statistic, choices = c("score","LR", "max"))
    if(any(exposure %in% names(coef(x)) == FALSE)){
        stop("exposure does not correspond to a coefficient in \'x\' \n")
    }

    # {{{ define all possible links newlinks
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
        ## directive <- res.find$directional
        ## restricted <- res.find$M.links
        ## link <- res.find$link

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
    
    # }}}
      
    # {{{ arguments LVM
    add.args <- setdiff(names(x$call), c("","x","data","control"))
    ls.LVMargs <- lapply(add.args, function(arg){x$call[[arg]]})
    names(ls.LVMargs) <- add.args
    if(is.null(data)){
        ls.LVMargs$data <- x$data$model.frame
    }else{
        ls.LVMargs$data <- data
    }
    ls.LVMargs$control <- x$control
    ls.LVMargs$control$trace <- FALSE
    # }}}
    
    # {{{ initialisation
    if(is.null(nStep)){
        nStep <- length(link)
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
    # }}}
    # {{{ loop over the models
    while(iStep <= nStep && NROW(iRestricted)>0 && cv==FALSE){

        if(trace >= 1){cat("Step ",iStep,":\n",sep="")}
        ls.LVMargs$start <- coef(iObject)

        if(statistic == "score"){
            # {{{ run modelsearch
            res.search <- modelsearch(iObject, link = iLink, silent = (trace <= 1))

            iN.link <- length(iLink)        
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
                                         ls.LVMargs = ls.LVMargs, method.p.adjust = method.p.adjust,
                                        display.warnings = display.warnings, trace = trace-1)
    
            # }}}
        }else if(statistic == "max"){
            # {{{ run modelsearchMax
           res.search <- modelsearchMax(iObject, restricted = iRestricted, link = iLink, directive = iDirective, alpha = alpha,
                                         ls.LVMargs = ls.LVMargs, method.p.adjust = method.p.adjust, method.max = method.max,
                                         conditional = conditional, mu.conditional = mu.conditional, iid.conditional = iid.conditional,
                                         export.iid = max(1,export.iid), trace = trace-1, ncpus = ncpus)
            # }}}
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
            ls.LVMargs$x <- addLink(iObject$model, var1 = iRestricted[index.rm,1], var2 = iRestricted[index.rm,2],
                                    covariance = iDirective[index.rm])
            ls.LVMargs$start <- coef(iObject)
                        
             
            suppressWarnings(
                res.search$best.model <- tryCatch(do.call(estimate, args = ls.LVMargs),
                                                  error = function(x){NA},
                                                  finally = function(x){x})
            )
            if(res.search$best.model$opt$convergence>0){
                ls.LVMargs$control$start <- NULL
                suppressWarnings(
                    res.search$best.model <- tryCatch(do.call(estimate, args = ls.LVMargs),
                                                      error = function(x){NA},
                                                      finally = function(x){x})
                )
            }

            iObject <- res.search$best.model
            iLink <- iLink[-index.rm]
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
                                   n.sim = n.sim, ncpus = ncpus, trace = trace)
            
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
# }}}




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


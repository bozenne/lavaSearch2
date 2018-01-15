### Lava_modelsearchMax.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 30 2017 (18:32) 
## Version: 
## last-updated: jan 15 2018 (11:38) 
##           By: Brice Ozenne
##     Update #: 613
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - modelsearchMax
#' @title Testing the Relevance of Additional Links Using the Max Statistic
#' @description Testing the Relevance of Additional Links Using the Max Statistic.
#' 
#' @name modelsearchMax
#'
#' @return an object of class lvmfit
#'
#' @seealso \code{link{modelsearch2}}
#' 
#' @keywords internal
#'


## * Function - modelsearchMax
## ' @rdname modelsearchMax
modelsearchMax <- function(x, restricted, link, directive, packages,
                           update.FCT, update.args, iid.FCT,                           
                           alpha, method.p.adjust, method.max, n.sim = 1e3, 
                           iid.previous = NULL, quantile.previous = NULL, 
                           export.iid, trace, ncpus, initCpus){

    ## WARNING: do not put link as NULL for data.table since it is used as an argument by the function
    convergence <- p.value <- statistic <- NULL ## [:for CRAN check] data.table

    ### ** initialisation
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    n.link <- NROW(restricted)
    nObs <- NROW(update.args$data)

    best.test <- -Inf
    best.model <- NULL    
    iid.link <- NULL
    convergence <- rep(NA,n.link)

    typeSD <- attr(iid.FCT, "typeSD")
    df <- attr(iid.FCT, "df")
    adjust.residuals <- attr(iid.FCT, "adjust.residuals")
    prepareScore <- df
    
    ### ** wraper
    warper <- function(iterI){ # iterI <- 2

        out <- list(dt = data.table(statistic = as.numeric(NA),
                                    df = as.numeric(NA),
                                    p.value = as.numeric(NA),
                                    adjusted.p.value = as.numeric(NA),
                                    convergence = as.numeric(NA),
                                    coefBeta = as.numeric(NA)),
                    iid = NULL)
        ## *** fit new model
        newfit <- update.FCT(x, args = update.args,
                             restricted = restricted[iterI,], directive = directive[iterI])
        out$dt[1, c("convergence") := newfit$opt$convergence]
        
        ## *** extract influence function        
        if(class(newfit) != "try-error"){ # test whether the model was estimated
            if(newfit$opt$convergence == 0){ # test whether lvmfit has correctly converged

                ## extract coefficient
                new.coef <- stats::coef(newfit)
                if(link[iterI] %in% names(new.coef) == FALSE){
                    stop("Coefficient ",link[iterI]," not found \n",
                         "Possible coefficients: ",paste0(names(new.coef), collapse = " "),"\n")
                }
                out$dt[1, "coefBeta" := new.coef[link[iterI]]]
                
                ## extract degree of freedom and standard error
                if(df){
                    dVcov2(newfit, return.score = TRUE) <- adjust.residuals
                    out$iid <- (attr(newfit$dVcov, "score") %*% attr(newfit$dVcov, "vcov.param")[,link[iterI],drop=FALSE])
                    
                    C <- matrix(0,
                                nrow = 1,
                                ncol = length(new.coef),
                                dimnames = list(NULL, names(new.coef)))
                    C[1,link[iterI]] <- 1
                    
                    e.df <- lTest(newfit, C = C, adjust.residuals = adjust.residuals,
                                  Ftest = FALSE)
                    
                    out$dt[1, "df" := e.df[1, "df"]]

                    if(typeSD == "information"){
                        sd.coef <- e.df[1, "std"]
                        out$iid <- out$iid * sd.coef / sd(out$iid, na.rm = TRUE)
                    }else{
                        sd.coef <- stats::sd(out$iid, na.rm = TRUE)
                    }
                    
                }else{

                    out$iid <- sqrt(nObs)*iid.FCT(newfit)[,link[iterI],drop=FALSE]
                    if(typeSD == "information"){
                        ## NOTE: it assumes that whenever df is FALSE then adjust.residuals is FALSE
                        sd.coef <- sqrt(stats::vcov(newfit)[link[iterI],link[iterI]])
                        out$iid <- out$iid * sd.coef / sd(out$iid, na.rm = TRUE)
                    }else{
                        sd.coef <- stats::sd(out$iid, na.rm = TRUE)
                    }
                        
                }
                ## compute test statistic
                out$dt[1, "statistic" := abs(.SD$coefBeta/sd.coef)] ## keep .SD for clarity
            }
        }
        return(out)
    }
    
### ** get influence function
    if(trace>0){
        cat("gather influence functions \n")
    }
            
    if(ncpus>1){

        FCTcombine <- function(res1,res2){
            res <- list(dt = rbind(res1$dt,res2$dt),
                        iid = cbind(res1$iid,res2$iid))
            return(res)
        }

        if(initCpus){
            cl <- parallel::makeCluster(ncpus)
            doParallel::registerDoParallel(cl)
        }
    
        if(trace > 0){
            parallel::clusterExport(cl, varlist = "trace")
        }

        vec.packages <- c("lavaSearch2", "data.table", packages)
        i <- NULL # [:for CRAN check] foreach
        res <- foreach::`%dopar%`(
                            foreach::foreach(i = 1:n.link, .packages =  vec.packages,
                                             # .export = c("ls.LVMargs"),
                                             .combine = FCTcombine),
                            {
                                if(trace){
                                    if(!exists("pb")){
                                        pb <- tcltk::tkProgressBar("modelsearchMax:", min=1, max=n.link)
                                    }
                                    tcltk::setTkProgressBar(pb, i)
                                }
                                return(warper(i))
                            })

        if(initCpus){
            parallel::stopCluster(cl)
        }
        
    }else{

        if(trace>0){
            requireNamespace("pbapply")
            resApply <- pbapply::pblapply(1:n.link, warper)
        }else{
            resApply <- lapply(1:n.link, warper)
        }
        res <- list(dt = data.table::rbindlist(lapply(resApply,"[[","dt")),
                    iid = do.call(cbind,lapply(resApply,"[[","iid")))
        
    }
    dt.test <- cbind(link = link, res$dt)    
    iid.link <- res$iid

    if(all(dt.test$convergence!=0)){
        stop("none of the extended model has converged \n",
             "the additional links may be misspecified \n")
    }
    
### ** p.value
    indexCV <- dt.test[, .I[convergence==0]]

    if(df){
        dt.test[indexCV, "p.value" := 2*(1-stats::pt(abs(.SD$statistic), df = df))] ## keep .SD for clarity
    }else{
        dt.test[indexCV, "p.value" := 2*(1-pnorm(abs(statistic)))]
    }

    ### ** adjust p.value
    if(method.p.adjust == "fastmax"){

        adj.tempo <- stats::p.adjust(dt.test$p.value, method = "bonferroni")
        if(any(adj.tempo<0.05)){
            dt.test[, p.value := as.numeric(p.value != min(p.value))]
            method.p.adjust <- "bonferroni"
        }else{
            method.p.adjust <- "max"
        }
        
    }
    
    if(method.p.adjust == "max"){
        nameN0 <- dt.test[indexCV, link]
        statisticN0 <- setNames(dt.test[convergence==0][["statistic"]],nameN0)
        if(df){
            dfN0 <- round(stats::median(stats::setNames(dt.test[convergence==0][["df"]],nameN0)))
        }else {
            dfN0 <- NULL
        }

        if(method.max=="integration"){
            resQmax <- calcDistMaxIntegral(statistic = statisticN0, iid = iid.link, df = dfN0,
                                           iid.previous = iid.previous, quantile.previous = quantile.previous, 
                                           alpha = alpha, ncpus = ncpus, initCpus = FALSE, trace = trace)
            
        }else{
            method.boot <- switch(method.max,
                                  "boot-naive" = "naive",
                                  "boot-residual" = "residual",
                                  "boot-wild" = "wild")
            
            resQmax <- calcDistMaxBootstrap(statistic = statisticN0, iid = iid.link, method = method.boot, n.sim = n.sim,
                                            iid.previous = iid.previous, quantile.previous = quantile.previous, 
                                            alpha = alpha, ncpus = ncpus, initCpus = FALSE, trace = trace)
        }

        dt.test[indexCV,c("corrected.level") := resQmax$correctedLevel]
        dt.test[indexCV,c("adjusted.p.value") := resQmax$p.adjust]
        dt.test[indexCV,c("quantile") := resQmax$z]
        Sigma <- resQmax$Sigma
        rownames(Sigma) <- dt.test[indexCV, link]
        colnames(Sigma) <- dt.test[indexCV, link]
        
        if(initCpus){
            parallel::stopCluster(cl)
        }
        
    }else{
        dt.test[indexCV,c("corrected.level") := NA]
        dt.test[dt.test$convergence==0,
                c("adjusted.p.value") := stats::p.adjust(p.value, method = method.p.adjust)]
        dt.test[indexCV,c("quantile") := as.numeric(NA)]
        Sigma <- NULL        
    }    

    ### ** export
    out <- list(dt.test = dt.test,
                iid = if(export.iid){iid.link}else{NULL},
                Sigma = Sigma)
    return(out)
}

## ----------------------------------------------------------------------
## Lava_modelsearchMax.R ends here

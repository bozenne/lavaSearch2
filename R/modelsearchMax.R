### Lava_modelsearchMax.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 30 2017 (18:32) 
## Version: 
## last-updated: okt  5 2017 (09:06) 
##           By: Brice Ozenne
##     Update #: 472
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - modelsearchMax
## ' @title Model searching using the max statistic
## ' @description Model searching using the max statistic to retain or not a link
## ' 
## ' @name modelsearchMax
## '
## ' @return an object of class lvmfit
## '
## ' @seealso \code{link{modelsearch2}}
## ' 
## ' @keywords internal
## '


## * Function - modelsearchMax
## ' @rdname modelsearchMax
modelsearchMax <- function(x, restricted, link, directive, packages,
                           update.FCT, update.args, iid.FCT,
                           method.p.adjust, method.max = "integration", n.sim = 1e3, alpha = 0.05,  
                           iid.previous = NULL, quantile.previous = NULL, 
                           export.iid = 1, trace = 1, ncpus = 1, initCpus = TRUE){

    ## ** initialisation
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    n.link <- NROW(restricted)
    nObs <- NROW(update.args$data)

    best.test <- -Inf
    best.model <- NULL    
    iid.link <- NULL
    convergence <- rep(NA,n.link)


    ## ** wraper
    warper <- function(iterI){ # iterI <- 2
        out <- list(dt = data.table(statistic = as.numeric(NA),
                                    p.value = as.numeric(NA),
                                    adjusted.p.value = as.numeric(NA),
                                    convergence = as.numeric(NA),
                                    coefBeta = as.numeric(NA)),
                    iid = NULL)
        ## *** fit new model
        newfit <- update.FCT(x, args = update.args, restricted = restricted[iterI,], directive = directive[iterI])
        out$dt[1, c("convergence") := newfit$opt$convergence]

        ## *** extract influence function        
        if(class(newfit) != "try-error"){ # test whether the model was estimated
            if(newfit$opt$convergence == 0){ # test whether lvmfit has correctly converged
                new.coef <- coef(newfit)
                if(link[iterI] %in% names(new.coef) == FALSE){
                    stop("Coefficient ",link[iterI]," not found \n",
                         "Possible coefficients: ",paste0(names(new.coef), collapse = " "),"\n")
                }
                new.coef <- new.coef[link[iterI]]
                #sd.coef <- sqrt(vcov(newfit)[link[iterI],link[iterI]])
                rdf <- df.residual(newfit, conservative = FALSE)

            out$dt[1, c("coefBeta") := new.coef]
            # test <- iidJack(newfit)

                out$iid <- sqrt(nObs)*iid.FCT(newfit)[,link[iterI],drop=FALSE]
                SeBeta <- sd(out$iid,na.rm = TRUE)
                out$dt[1, c("statistic") := abs(.SD$coefBeta/SeBeta)]               
            }
        }
        return(out)
    }
    
    ## ** get influence function
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
            doSNOW::registerDoSNOW(cl)
        }
    
        if(trace > 0){
            pb <- utils::txtProgressBar(max = n.link, style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
        }else{
            opts <- NULL
        }

        vec.packages <- c("lavaSearch2", "data.table", packages)
        i <- NULL # for CRAN check
        res <- foreach::`%dopar%`(
                            foreach::foreach(i = 1:n.link, .packages =  vec.packages,
                                             # .export = c("ls.LVMargs"),
                                             .combine = FCTcombine,
                                             .options.snow = opts),
                            {
                                return(warper(i))
                            })
   
        if(trace > 0){  close(pb) }
        
    }else{
        
        if(trace>0){
            requireNamespace("pbapply")
            coefJack <- pbapply::pblapply(1:n.link, warper)
        }else{
            coefJack <- lapply(1:n.link, warper)
        }
        res <- list(dt = data.table::rbindlist(lapply(coefJack,"[[","dt")),
                    iid = do.call(cbind,lapply(coefJack,"[[","iid")))
        
    }
    dt.test <- cbind(link = link, res$dt)    
    iid.link <- res$iid
    
### ** p.value
    df.model <- df.residual(x, conservative = TRUE)
    indexCV <- dt.test[, .I[.SD$convergence==0]]

    if(is.null(df.model)){
        dt.test[indexCV, c("p.value") := 2*(1-pnorm(abs(.SD$statistic)))]
    }else{
        dt.test[indexCV, c("p.value") := 2*(1-pt(abs(.SD$statistic), df = df.model))]
    }

    ### ** adjust p.value
    if(method.p.adjust == "max"){
        nameN0 <- dt.test[indexCV, .SD$link]
        statisticN0 <- setNames(dt.test[convergence==0][["statistic"]],nameN0)

        if(method.max=="integration"){
            args(calcDistMaxIntegral)
            resQmax <- calcDistMaxIntegral(statistic = statisticN0, iid = iid.link, df = df.model,
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
        rownames(Sigma) <- dt.test[indexCV,.SD$link]
        colnames(Sigma) <- dt.test[indexCV,.SD$link]
        
        if(initCpus){
            parallel::stopCluster(cl)
        }
        
    }else{
        dt.test[indexCV,c("corrected.level") := NA]
        dt.test[dt.test$convergence==0, c("adjusted.p.value") := p.adjust(.SD$p.value, method = method.p.adjust)]
        dt.test[indexCV,c("quantile") := as.numeric(NA)]
        Sigma <- NULL        
    }    
    
## ** export
    out <- list(dt.test = dt.test,
                iid = if(export.iid){iid.link}else{NULL},
                Sigma = Sigma)
    return(out)
}

## ----------------------------------------------------------------------
## Lava_modelsearchMax.R ends here

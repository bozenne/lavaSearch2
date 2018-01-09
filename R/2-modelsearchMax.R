### Lava_modelsearchMax.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 30 2017 (18:32) 
## Version: 
## last-updated: jan  8 2018 (15:37) 
##           By: Brice Ozenne
##     Update #: 533
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - modelsearchMax
#' @title Model searching using the max statistic
#' @description Model searching using the max statistic to retain or not a link
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
                           robust = FALSE, df = TRUE, adjust.residuals = TRUE,
                           method.p.adjust, method.max = "integration", n.sim = 1e3, alpha = 0.05,  
                           iid.previous = NULL, quantile.previous = NULL, 
                           export.iid = 1, trace = 1, ncpus = 1, initCpus = TRUE){

### ** initialisation
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    n.link <- NROW(restricted)
    nObs <- NROW(update.args$data)

    best.test <- -Inf
    best.model <- NULL    
    iid.link <- NULL
    convergence <- rep(NA,n.link)

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
                new.coef <- coef(newfit)
                if(link[iterI] %in% names(new.coef) == FALSE){
                    stop("Coefficient ",link[iterI]," not found \n",
                         "Possible coefficients: ",paste0(names(new.coef), collapse = " "),"\n")
                }
                out$dt[1, "coefBeta" := new.coef[link[iterI]]]

                ## extract degree of freedom and standard error
                if(df){
                    C <- matrix(0,
                                nrow = 1,
                                ncol = length(new.coef),
                                dimnames = list(NULL, names(new.coef)))
                    C[1,link[iterI]] <- 1
                    e.df <- lTest(newfit, C = C, adjust.residuals = adjust.residuals, Ftest = FALSE)
                    out$dt[1, "df" := e.df[1, "df"]]

                    sd.coef <- e.df[1, "std"]
                }else{
                    if(robust == FALSE){
                        sd.coef <- sqrt(vcov(newfit)[link[iterI],link[iterI]])
                    }
                }
                
                ## extract iid
                out$iid <- sqrt(nObs)*iid.FCT(newfit, adjust.residuals = adjust.residuals)[,link[iterI],drop=FALSE]
                if(robust == FALSE){                    
                    out$iid <- out$iid * sd.coef / sqrt(mean(out$iid^2, na.rm = TRUE))
                }else{
                    sd.coef <- sd(out$iid, na.rm = TRUE)
                }

                ## compute test statistic
                out$dt[1, "statistic" := abs(.SD$coefBeta/sd.coef)]
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
    indexCV <- dt.test[, .I[.SD$convergence==0]]

    if(df){
        dt.test[indexCV, "p.value" := 2*(1-pnorm(abs(.SD$statistic)))]
    }else{
        dt.test[indexCV, "p.value" := 2*(1-pt(abs(.SD$statistic), df = df))]
    }

    ### ** adjust p.value
    if(method.p.adjust == "max"){
        nameN0 <- dt.test[indexCV, .SD$link]
        statisticN0 <- setNames(dt.test[convergence==0][["statistic"]],nameN0)
        if(df){
            dfN0 <- round(median(setNames(dt.test[convergence==0][["df"]],nameN0)))
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
    
### ** export
    out <- list(dt.test = dt.test,
                iid = if(export.iid){iid.link}else{NULL},
                Sigma = Sigma)
    return(out)
}

## ----------------------------------------------------------------------
## Lava_modelsearchMax.R ends here

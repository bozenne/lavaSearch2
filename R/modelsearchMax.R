### Lava_modelsearchMax.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 30 2017 (18:32) 
## Version: 
## last-updated: jun 23 2017 (19:58) 
##           By: Brice Ozenne
##     Update #: 265
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

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

#' @rdname modelsearchMax
modelsearchMax <- function (x, data = NULL, restricted, link, directive, method.max = "integration", method.iid = "iidJack", alpha = 0.05,
                            ls.LVMargs = NULL, method.p.adjust = "bonferroni", n.sim = 1e4,
                            conditional = NULL, mu.conditional = NULL, iid.conditional = NULL, 
                            export.iid = 1, trace = 1, ncpus = 1, initCpus = TRUE){

    `p.value (Wald)` <- `statistic (Wald)` <- statistic <- p.value <- adjusted.p.value <- `adjusted.p.value (Wald)` <- coefBeta <- NULL # for CRAN check

                                        # {{{ initialisation
    if(is.null(ncpus)){ ncpus <- parallel::detectCores()}
    n.link <- length(link)
    nObs <- x$data$n
    dt.test <- data.table("link" = link,
                          "statistic" = as.numeric(rep(NA,n.link)),
                          "p.value" = as.numeric(rep(NA,n.link)),
                          "adjusted.p.value" = as.numeric(rep(NA,n.link)),
                          "statistic (Wald)" = as.numeric(rep(NA,n.link)),
                          "p.value (Wald)" = as.numeric(rep(NA,n.link)),
                          "adjusted.p.value (Wald)" = as.numeric(rep(NA,n.link))
                          )

    

    best.test <- -Inf
    best.model <- NULL    
    iid.link <- NULL
    cv <- rep(NA,n.link)
    # }}}

   
                                        # {{{ wraper                                 
    warper <- function(iterI){ # iterI <- 2
        out <- list(dt = data.table(statistic = as.numeric(NA),
                                    p.value = as.numeric(NA),
                                    adjusted.p.value = as.numeric(NA),
                                    `statistic (Wald)` = as.numeric(NA),
                                    `p.value (Wald)` = as.numeric(NA),
                                    `adjusted.p.value (Wald)` = as.numeric(NA),
                                    cv = as.numeric(NA),
                                    coefBeta = as.numeric(NA)),
                    iid = NULL)

        # {{{ fit new model
        ls.LVMargs$x <- addLink(x$model, var1 = restricted[iterI,1], var2 = restricted[iterI,2],
                                covariance = (1-directive[iterI]))

        suppressWarnings(
            newfit <- tryCatch(do.call(estimate, args = ls.LVMargs),
                               error = function(x){NA},
                               finally = function(x){x})
        )
        if(newfit$opt$convergence>0){
            ls.LVMargs$control$start <- NULL
            suppressWarnings(
                newfit <- tryCatch(do.call(estimate, args = ls.LVMargs),
                                   error = function(x){NA},
                                   finally = function(x){x})
            )
        }
        out$dt[1,cv := newfit$opt$convergence]
        # }}}

        # {{{ extract influence function
        if("lvmfit" %in% class(newfit)){ # test lvmfit is not an error
            if(newfit$opt$convergence == 0){ # test whether lvmfit has correctly converged
               # if(restricted[iterI,2] %in% endogenous(x) && restricted[iterI,1] %in% endogenous(x)){
               #     nameNewlink <- restricted[iterI,2]
               # }else{
                 nameNewlink <- link[iterI]
                                        # }
                out$dt[1, coefBeta := coef(newfit)[nameNewlink]]
                out$iid <- sqrt(nObs)*do.call(method.iid, args = list(newfit, cpus = 1))[,nameNewlink,drop=FALSE]
                SeBeta <- sd(out$iid,na.rm=TRUE) # note n/n-1 vs. sqrt(vcov(newfit)[iLink[iterI],iLink[iterI]])
                out$dt[1,`statistic` := abs(coefBeta/SeBeta)]

                #IF.beta <- sqrt(nObs)*iid(newfit)[,nameNewlink,drop=FALSE]
                #SeBeta <- sd(IF.beta)

                out$dt[1,`statistic (Wald)` := summary(newfit)$coef[nameNewlink,1]/summary(newfit)$coef[nameNewlink,2]]
                out$dt[1,`p.value (Wald)` := summary(newfit)$coef[nameNewlink,4]]                
            }
        }
                                        # }}}
        return(out)
    }
                                        # }}}

                                        #  sqrt(nObs)*
                                        # {{{ parallel computations
    if(initCpus){
        cl <- parallel::makeCluster(ncpus)
        doSNOW::registerDoSNOW(cl)
    }
    
    if(trace > 0){
        cat("gather influence functions \n")
        pb <- utils::txtProgressBar(max = n.link, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
    }else{
        opts <- NULL
    }

    FCTcombine <- function(res1,res2){
        res <- list(dt = rbind(res1$dt,res2$dt),
                    iid = cbind(res1$iid,res2$iid))
        return(res)
    }
    vec.packages <- c("lava", "lavaSearch2", "data.table")
    i <- NULL
    res <- foreach::`%dopar%`(
                        foreach::foreach(i = 1:n.link, .packages =  vec.packages,
                                         .export = c("ls.LVMargs"),
                                         .combine = FCTcombine,
                                         .options.snow = opts),
                        {
                            return(warper(i))
                        })
   
    dt.test <- cbind(link = link, res$dt)
    iid.link <- res$iid
    
    if(trace > 0){  close(pb) }                                           
    # }}}
    
    ## if(na.omit){
    ##     index.NNA <- which(rowSums(!is.na(dt.test))>0)
    ##     dt.test <- dt.test[index.NNA]
    ## }

    
                                        # {{{ adjust p.value
    dt.test[cv==0, `adjusted.p.value (Wald)` := p.adjust(`p.value (Wald)`, method = method.p.adjust)]
    dt.test[cv==0, p.value := 2*(1-pnorm(abs(statistic)))]

    iid.all <- cbind(iid.conditional,iid.link)
    index.NNA <- which(rowSums(is.na(iid.all))==0)
    iid.all <- iid.all[index.NNA,,drop=FALSE]
    mu.all <- c(mu.conditional,dt.test[cv==0,coefBeta])
    mu.conditional.all <- c(conditional,rep(0,n.link))
    resQmax <- calcDistMax(dt.test[cv==0,statistic], iid = iid.all, mu = mu.all,
                           df = degreeOfFreedom(x, conservative = TRUE)-1, alpha = alpha,
                           conditional = mu.conditional.all, method = method.max,
                           n.sim = n.sim, ncpus = ncpus, initCpus = FALSE, trace = trace)
   
    dt.test[cv==0,adjusted.p.value := resQmax$p.adjust]
    z <- resQmax$z

    
    if(initCpus){
        parallel::stopCluster(cl)
    } 
    # }}}
    
    #### export
    name.max <- dt.test[which.max(abs(dt.test[["statistic"]])),link]
    
    out <- list(dt.test = dt.test,
                iid.all = if(export.iid>1){iid.link}else{NULL},
                iid.link = if(export.iid>0){iid.link[,name.max,drop=FALSE]}else{NULL},
                z = z)
    return(out)
}

#----------------------------------------------------------------------
### Lava_modelsearchMax.R ends here

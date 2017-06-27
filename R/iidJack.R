### iidJack.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 23 2017 (09:15) 
## Version: 
## last-updated: jun 27 2017 (11:47) 
##           By: Brice Ozenne
##     Update #: 67
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ Documentation
#' @title Jacknife i.i.d. decomposition (influence function) from model object
#' @description Extract i.i.d. decomposition (influence function) from model object
#'
#' @name iidJack
#' 
#' @param x model object.
#' @param data dataset used to perform the jacknife.
#' @param ncpus number of cpus available for parallel computation.
#' @param initCpus should the parallel computation be initialized?
#' @param trace should a progress bar be used to trace the execution of the function
#' @param ... additional arguments.
#' @examples
#' n <- 20
#'
#' #### glm ####
#' m <- lvm(y~x+z)
#' distribution(m, ~y+z) <- binomial.lvm("logit")
#' d <- sim(m,n)
#' g <- glm(y~x+z,data=d,family="binomial")
#' iidJack(g)
#'
#' #### Cox model ####
#' library(survival)
#' library(riskRegression)
#' data(Melanoma, package = "riskRegression")
#' m <- coxph(Surv(time,status==1)~ici+age, data = Melanoma, x = TRUE, y = TRUE)
#' iid1 <- iidJack(m)
#' iid2 <- iidCox(m)$IFbeta
#' 
#' apply(iid1,2,sd)
#' apply(iid2,2,sd)
#'
#' #### LVM ####
#' set.seed(10)
#' m <- lvm(Y~X1+X2+X3)
#' d <- sim(m, n)
#' e <- estimate(m, data = d)
#' iid1 <- iidJack(e)
#' iid2 <- iid(e)
#' attr(iid2, "bread") <- NULL
#'
#' apply(iid1,2,sd)
#' apply(iid2,2,sd)
#' quantile(iid2 - iid1)
#' 
#' @export
iidJack <- function(x,...) UseMethod("iidJack")
# }}}

#' @rdname iidJack
#' @export
iidJack.default <- function(x,data=NULL,ncpus=1,initCpus,trace=TRUE,...) {
    
    estimate.lvm <- lava_estimate.lvm
    
    # {{{ extract data
    if(is.null(data)){
        if("lvmfit" %in% class(x)){
            data <- model.frame(x, all = TRUE)
        }else{
            data <- eval(x$call$data)
        }
    }
    n.obs <- NROW(data)
    coef.x <- coef(x)
    # }}}

    # {{{ warper
    warper <- function(i){ # i <- 1
        xnew <- update(x, data = data[-i,])
        return(c(coef.x-coef(xnew)))
        ## return(c(coef(xnew)-coef.x,
        ##          mu = predict(xnew, data = data[i,,drop=FALSE]))
        ##        )
    }
    # }}}
    
    # {{{ parallel computations
    if(ncpus>1){
        if(initCpus){
            cl <- parallel::makeCluster(ncpus)
            doSNOW::registerDoSNOW(cl)
        }
    
        if(trace > 0){
            cat("jacknife \n")
            pb <- utils::txtProgressBar(max = n.obs, style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
        }else{
            opts <- NULL
        }

        vec.packages <- c("lava")
        toExport <- sapply(as.list(x$call),as.character)
        i <- NULL
        coefJack <- foreach::`%dopar%`(
                                 foreach::foreach(i = 1:n.obs, .packages =  vec.packages,
                                                  .export = toExport,
                                                  .combine = "rbind",
                                                  .options.snow = opts),{
                                                      warper(i)
                                                  })
    
    
        if(initCpus){
            parallel::stopCluster(cl)
        }

        if(trace > 0){ close(pb) }                                           
        rownames(coefJack) <- 1:n.obs
    }else{
        coefJack <- lapply(1:n.obs, warper)
        coefJack <- do.call(rbind, coefJack)
    }
    # }}}

    return(coefJack)
}
    
#----------------------------------------------------------------------
### iidJack.R ends here

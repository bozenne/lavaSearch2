### iidJack.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 23 2017 (09:15) 
## Version: 
## last-updated: aug 28 2017 (09:28) 
##           By: Brice Ozenne
##     Update #: 146
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
#' @param grouping variable defining cluster of observations that will be simultaneously removed by the jackknife.
#' @param ncpus number of cpus available for parallel computation.
#' @param initCpus should the parallel computation be initialized?
#' @param trace should a progress bar be used to trace the execution of the function
#' @param ... additional arguments.
#' @examples
#' n <- 20
#'
#' #### glm ####
#' set.seed(10)
#' m <- lvm(y~x+z)
#' distribution(m, ~y+z) <- binomial.lvm("logit")
#' d <- sim(m,n)
#' g <- glm(y~x+z,data=d,family="binomial")
#' iid1 <- iidJack(g)
#' iid2 <- iid(g)
#' quantile(iid1-iid2)
#' vcov(g)
#' colSums(iid2^2)
#' colSums(iid1^2)
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
#'
#' print(iid2)
#' 
#' apply(iid2,2,sd)
#'
#' #### LVM ####
#' set.seed(10)
#'
#' mSim <- lvm(c(Y1,Y2,Y3,Y4,Y5) ~ 1*eta)
#' latent(mSim) <- ~eta
#' categorical(mSim, K=2) <- ~G
#' dW <- as.data.table(sim(mSim, n, latent = FALSE))
#' dW[,Id := as.character(1:.N)]
#' dL <- melt(dW, id.vars = c("G","Id"), variable.name = "time", value.name = "Y")
#' dL[,time := gsub("Y","",time)]
#'
#' m1 <- lvm(c(Y1,Y2,Y3,Y4,Y5) ~ 1*eta)
#' latent(m1) <- ~eta
#' regression(m1) <- eta ~ G
#' e <- estimate(m1, data = dW)
#' iid1 <- iidJack(e)
#' iid2 <- iid(e)
#' attr(iid2, "bread") <- NULL
#'
#' apply(iid1,2,sd)
#' apply(iid2,2,sd)
#' quantile(iid2 - iid1)
#'
#' library(nlme)
#' e2 <- lme(Y~G+time, random = ~1|Id, weights = varIdent(form =~ 1|Id), data = dL)
#' e2 <- lme(Y~G, random = ~1|Id, data = dL)
#' iid3 <- iidJack(e2)
#' apply(iid3,2,sd)
#'
#' @export
iidJack <- function(x,...) UseMethod("iidJack")
# }}}

#' @rdname iidJack
#' @export
iidJack.default <- function(x,data=NULL,grouping=NULL,ncpus=1,initCpus=TRUE,trace=TRUE,...) {
    
    estimate.lvm <- lava_estimate.lvm
    
    # {{{ extract data
    if(is.null(data)){
        if("lvmfit" %in% class(x)){
            data <- model.frame(x, all = TRUE)
        }else{
            data <- eval(x$call$data)
        }
    }
    if(is.data.table(data)){
        data <- copy(data)
    }else{
        data <- as.data.table(data)
    }
    n.obs <- NROW(data)
    if(any(class(x) %in% "lme")){
        getCoef <- nlme::fixef
    }else{
        getCoef <- coef
    }
    coef.x <- getCoef(x)
    names.coef <- names(coef.x)
    n.coef <- length(coef.x)
    # }}}

    # {{{ define the grouping level for the data
    if(is.null(grouping)){
        if(any(class(x)%in%c("lme","gls","nlme"))){
            data[, c("XXXgroupingXXX") := as.vector(apply(x$groups,2,interaction))]
        }else{
            data[, c("XXXgroupingXXX") := 1:NROW(data)]
        }
        grouping <- "XXXgroupingXXX"        
    }else{
        if(length(grouping)>1){
            stop("grouping must refer to only one variable \n")
        }
        if(grouping %in% names(data) == FALSE){
            stop("variable defined in grouping not found in data \n")
        }
    }
    data[, c(grouping) := as.character(.SD$grouping)]
    Ugrouping <- unique(data[[grouping]])
    n.group <- length(Ugrouping)
    # }}}
    
    # {{{ warper
    warper <- function(i){ # i <- 1
        xnew <- try(update(x, data = data[data[[grouping]]!=i,]), silent = TRUE)
        if(class(xnew)!="try-error"){
            return(getCoef(xnew))
        }else{
            return(rep(NA, n.coef))
        }
        ## return(c(coef(xnew)-coef.x,
        ##          mu = predict(xnew, data = data[i,,drop=FALSE]))
        ##        )
    }
    # }}}
    
    # {{{ parallel computations: get jackknife coef
   ## if(ncpus>1){
        if(initCpus){
            cl <- parallel::makeCluster(ncpus)
            doSNOW::registerDoSNOW(cl)
        }
    
        if(trace > 0){
            cat("jacknife \n")
            pb <- utils::txtProgressBar(max = n.group, style = 3)
            progress <- function(n) setTxtProgressBar(pb, n)
            opts <- list(progress = progress)
        }else{
            opts <- NULL
        }


    estimator <- as.character(x$call[[1]]) 

    vec.packages <- c("lava")
    possiblePackage <- gsub("package:","",utils::getAnywhere(estimator)$where[1])
    existingPackage <- as.character(utils::installed.packages()[,"Package"])

    ls.call <- as.list(x$call)
    test.length <- which(unlist(lapply(ls.call, length))==1)
    test.class <- which(unlist(lapply(ls.call, function(cc){
                                 (class(c) %in% c("numeric","character","logical")) == FALSE
    })))
    test.class <- which(unlist(lapply(ls.call, class)) %in% c("numeric","character","logical") == FALSE)
    
    indexExport <- intersect(test.class,test.length)
    toExport <- sapply(ls.call[indexExport], as.character)
    
    if(possiblePackage %in% existingPackage){
        vec.packages <- c(vec.packages,possiblePackage)
    }
    if(length(x$call$data)==1){
        toExport <- c(toExport,as.character(x$call$data))
    }
    if(length(x$call$formula)==1){
        toExport <- c(toExport,as.character(x$call$formula))
    }
    if(length(x$call$fixed)==1){
        toExport <- c(toExport,as.character(x$call$fixed))        
    }

    #sapply(as.list(x$call),as.character)
    i <- NULL # for CRAN check
    coefJack <- foreach::`%dopar%`(
                             foreach::foreach(i = Ugrouping, .packages =  vec.packages,
                                              .export = toExport,
                                              .combine = "rbind",
                                              .options.snow = opts),{
                                                  warper(i)
                                              })
    
    if(initCpus){
        parallel::stopCluster(cl)
    }
    if(trace > 0){ close(pb) }                                           
    rownames(coefJack) <- 1:n.group
       ## }else{
       ##     coefJack <- lapply(Ugrouping, warper)
       ##     coefJack <- do.call(rbind, coefJack)
       ## }
    # }}}

    # {{{ post treatment: from jackknife coef to IF
    # defined as (n-1)*(coef-coef(-i))
    # division by n to match output of lava, i.e. IF/n
    iidJack <- -(n.group-1)/n.group*sweep(coefJack, MARGIN = 2, STATS = coef.x, FUN = "-")
    colnames(iidJack) <- names.coef
    # }}}

    return(iidJack)
}
    
#----------------------------------------------------------------------
### iidJack.R ends here

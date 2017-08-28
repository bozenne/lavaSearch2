### compareIID.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jul 18 2017 (17:15) 
## Version: 
## last-updated: jul 18 2017 (17:59) 
##           By: Brice Ozenne
##     Update #: 25
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @param model a lvm model .
#' @param data a dataset.
#' @param condition the variable in data indicating the different conditions to be compared.
#' @param model.init the model used to initialize the coefficient.
#' @param control control argument of estimate.
#' @param ... additional argument to be passed to estimate.
#' 
#' @examples
#' 
#' library(lava)
#' 
#' m <- lvm(Y~X1+X2+X3)
#' m0 <- lvm(Y ~ 1*X1+1*X2+1*X3)
#' m1 <- lvm(Y ~ 2*X1+1*X2+1*X3)
#' n <- 1e2  
#'  
#' ## under H0
#' set.seed(1)
#' d1 <- sim(m0, 1e2)
#' d2 <- sim(m0, 1e2)
#' d <- rbind(cbind(d1,cond = 1), cbind(d2,cond = 2))
#' compareCondition(m, data = d, condition = "cond")
#' 
#' ## under H1
#' set.seed(1)
#' d1 <- sim(m0, 1e2)
#' d2 <- sim(m1, 1e2)
#' d <- rbind(cbind(d1,cond = 1), cbind(d2,cond = 2))
#' compareCondition(m, data = d, condition = "cond")
#' 
compareCondition <- function(model, data, condition, coef = NULL,
                             model.init = NULL, control = list(), trace = TRUE, ...){

    ## normalize data

    data <- copy(as.data.table(data))
    ## find all conditions
    Ucondition <- unique(data[[condition]])
    n.condition <- length(Ucondition)

    ## estimate all models
    if(trace){cat("Estimate all models: ")}
    ls.e <- list()
    ls.beta <- list()
    ls.iid <- list()
    for(iCondition in 1:n.condition){ # iCondition <- 1
        indexCond <- which(data[[condition]]==Ucondition[iCondition])
        iData <- data[indexCond]

        if(!is.null(model.init) && "start" %in% names(control) == FALSE){
            control$start <- coef(estimate(model, data = iData, control = control))            
        }
        ls.e[[iCondition]] <- estimate(model, data = iData,
                                       control = control, ...)
        ls.beta[[iCondition]] <- data.table(Ucondition[iCondition],
                                            ls.e[[iCondition]]$data$n,
                                            coef(ls.e[[iCondition]]),
                                            confint(ls.e[[iCondition]]))
        names(ls.beta[[iCondition]]) <- c("condition","n","Estimate","ci.inf","ci.sup")
    }
    names(ls.beta) <- Ucondition
    names(ls.e) <- Ucondition
    if(trace){cat("done\n")}
    dt.beta <- rbindlist(ls.beta)

    ## pairwise comparison
    grid <- expand.grid(1:n.condition,1:n.condition)
    # remove same condition and symmetric
    grid <- grid[grid[,1]>grid[,2],]
    grid <- cbind(Ucondition[grid[,1]],Ucondition[grid[,2]])
    n.grid <- NROW(grid)

    dt.res <- NULL
    qMax <- rep(NA, n.grid)

    if(trace){cat("Estimate the distribution of the max statistic: ")}
    for(iGrid in 1:n.grid){ # iGrid <- 1
        C1 <- grid[iGrid,1]
        C2 <- grid[iGrid,2]
        df.C12 <- ls.e[[C1]]$data$n+ls.e[[C2]]$data$n
        
        ## coef to keep
        if(is.null(coef)){            
            iCoef <- names(coef(ls.e[[C1]]))
        }else{
            iCoef <- coef
        }
        n.coef <- length(iCoef)

        ## iid
        iid.C1 <- iid(ls.e[[C1]],
                      data = model.frame(ls.e[[C1]], data = data[condition %in% c(C1,C2)]))[,iCoef,drop=FALSE]
        attr(iid.C1, "bread") <- NULL
        iid.C2 <- iid(ls.e[[C2]],
                      data = model.frame(ls.e[[C1]], data = data[condition %in% c(C1,C2)]))[,iCoef,drop=FALSE]
        attr(iid.C2, "bread") <- NULL

        ## statistic test
        s.C1 <- summary(ls.e[[C1]])$coef
        s.C2 <- summary(ls.e[[C2]])$coef

        diff.coef <- s.C1[iCoef,"Estimate"] - s.C2[iCoef,"Estimate"]
        diff.t <- diff.coef/(s.C1[iCoef,"Std. Error"] + s.C2[iCoef,"Std. Error"])

        ## find the quantile of the max statistic
        resQmax <- calcDistMax(statistic = diff.t,
                               iid = iid.C1 - iid.C2,
                               mu = diff.coef,
                               conditional = rep(0, n.coef),
                               df = df.C12,
                               method = "integration",
                               alpha = 0.05,
                               ncpus = 1,
                               initCpus = TRUE,
                               n.sim = NA,
                               trace = FALSE,
                               n.repMax = 100)

        ## output results
        iDt.res <- data.table(coef = iCoef,
                             condition1 = C1,
                             value1 = s.C1[iCoef,"Estimate"],
                             condition2 = C2,
                             value2 = s.C2[iCoef,"Estimate"],
                             diff = diff.coef,
                             p.nadjust = 2*(1-pt(abs(diff.t), df = df.C12))
                             )
        iDt.res[,p.holm :=  p.adjust(p.nadjust, method = "holm")]
        iDt.res[,p.max := resQmax$p.adjust]

        dt.res <- rbind(dt.res, iDt.res)
        qMax[iGrid] <- resQmax$z
    }
    if(trace){cat("done\n")}

    ## export
    return(list(test = dt.res, quantile.max = qMax))
}

#----------------------------------------------------------------------
### compareIID.R ends here

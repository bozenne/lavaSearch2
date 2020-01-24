### compare2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 30 2018 (14:33) 
## Version: 
## Last-Updated: jan 24 2020 (17:16) 
##           By: Brice Ozenne
##     Update #: 666
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation - compare2
#' @title Test Linear Hypotheses with small sample correction
#' @description Test Linear Hypotheses using a multivariate Wald statistic from \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} objects.
#' Similar to \code{lava::compare} but with small sample correction (if any).
#' @name compare2
#'
#' @param object a \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} object.
#' @param df [logical] should the degree of freedoms of the Wald statistic be computed using the Satterthwaite correction?
#' Otherwise the degree of freedoms are set to \code{Inf}, i.e. a normal distribution is used instead of a Student's t distribution when computing the p-values.
#' @param ssc [logical] should the standard errors of the coefficients be corrected for small sample bias? Argument passed to \code{sCorrect}.
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param par [vector of characters] expression defining the linear hypotheses to be tested.
#' See the examples section. 
#' @param contrast [matrix] a contrast matrix defining the left hand side of the linear hypotheses to be tested.
#' @param robust [logical] should the robust standard errors be used instead of the model based standard errors?
#' @param rhs [vector] the right hand side of the linear hypotheses to be tested.
#' @param as.lava [logical] should the output be similar to the one return by \code{lava::compare}?
#' @param F.test [logical] should a joint test be performed?
#' @param conf.level [numeric 0-1] the confidence level of the confidence interval.
#' @param ...  [internal] only used by the generic method.
#'
#' @details The \code{par} argument or the arguments \code{contrast} and \code{null} (or equivalenty \code{rhs})
#' specify the set of linear hypotheses to be tested. They can be written:
#' \deqn{
#'   contrast * \theta = null
#' }
#' where \eqn{\theta} is the vector of the model coefficients. \cr
#' The \code{par} argument must contain expression(s) involving the model coefficients.
#' For example \code{"beta = 0"} or \code{c("-5*beta + alpha = 3","-alpha")} are valid expressions if alpha and beta belong to the set of model coefficients.
#' A contrast matrix and the right hand side will be generated inside the function. \cr
#' 
#' When directly specified, the contrast matrix must contain as many columns as there are coefficients in the model (mean and variance coefficients).
#' Each hypothesis correspond to a row in the contrast matrix. \cr
#'
#' The null vector should contain as many elements as there are row in the contrast matrix. \cr
#' 
#' Argument rhs and null are equivalent.
#' This redondance enable compatibility between \code{lava::compare}, \code{compare2}, \code{multcomp::glht}, and \code{glht2}.
#'
#' @seealso \code{\link{createContrast}} to create contrast matrices. \cr
#' \code{\link{sCorrect}} to pre-compute quantities for the small sample correction.
#' 
#' @return If \code{as.lava=TRUE} an object of class \code{htest}.
#' Otherwise a \code{data.frame} object.

## * example - compare2
#' @examples
#' #### simulate data ####
#' set.seed(10)
#' mSim <- lvm(Y~0.1*X1+0.2*X2)
#' categorical(mSim, labels = c("a","b","c")) <- ~X1
#' transform(mSim, Id~Y) <- function(x){1:NROW(x)}
#' df.data <- lava::sim(mSim, 1e2)
#'
#' #### with lm ####
#' ## direct use of compare2
#' e.lm <- lm(Y~X1+X2, data = df.data)
#' anova(e.lm)
#' compare2(e.lm, par = c("X1b=0","X1c=0"))
#' 
#' ## or first compute the derivative of the information matrix
#' sCorrect(e.lm) <- TRUE
#' 
#' ## and define the contrast matrix
#' C <- createContrast(e.lm, linfct = c("X1b=0","X1c=0"), add.variance = TRUE)
#'
#' ## run compare2
#' compare2(e.lm, linfct = C$contrast, rhs = C$null)
#' compare2(e.lm, linfct = C$contrast, rhs = C$null, robust = TRUE)
#' 
#' #### with gls ####
#' library(nlme)
#' e.gls <- gls(Y~X1+X2, data = df.data, method = "ML")
#'
#' ## first compute the derivative of the information matrix
#' sCorrect(e.gls, cluster = 1:NROW(df.data)) <- TRUE
#' 
#' compare2(e.gls, par = c("5*X1b+2*X2 = 0","(Intercept) = 0"))
#' 
#' #### with lvm ####
#' m <- lvm(Y~X1+X2)
#' e.lvm <- estimate(m, df.data)
#' 
#' compare2(e.lvm, par = c("-Y","Y~X1b+Y~X1c"))
#' compare2(e.lvm, par = c("-Y","Y~X1b+Y~X1c"), robust = TRUE)
#' @concept small sample inference
#' @export
`compare2` <-
    function(object, linfct, rhs, robust, cluster,
             ssc, df,
             as.lava, F.test, conf.level) UseMethod("compare2")

## * compare2.lm
#' @rdname compare2
#' @export
compare2.lm <- function(object, linfct = NULL, rhs = NULL,
                        robust = FALSE, cluster = NULL, 
                        ssc = lava.options()$ssc, df = lava.options()$df,
                        as.lava = TRUE, F.test = TRUE, conf.level = 0.95){

    if(is.null(linfct)){ ## necessary for lava::gof to work
        return(lava::compare(object))
    }
    
    if(is.null(object$sCorrect) || !identical(object$sCorrect$ssc$type, ssc) || !identical(object$sCorrect$df, df)){
        object <- sCorrect(object, ssc = ssc, df = df)
    }

    if(!is.logical(robust)){ 
        stop("Argument \'robust\' should be TRUE or FALSE \n")
    }

    if(robust){
        factor.dRvcov <- lava.options()$factor.dRvcov

        if(!is.null(cluster) && !inherits(cluster, "function")){
            
            if(length(cluster)==1){
                ## reconstruct cluster variable
                data <- object$sCorrect$data
                
                if(cluster %in%  names(data) == FALSE){
                    stop("Could not find variable ",cluster," (argument \'cluster\') in argument \'data\' \n")
                }else{
                    cluster <- data[[cluster]]
                }            
                
            }else if(NROW(data)!=length(cluster)){
                stop("length of argument \'cluster\' does not match number of rows of the dataset \n")
            }
            n.cluster <- length(unique(cluster))
        }else{
            n.cluster <- stats::nobs(object)
            cluster <- NULL
        }
    }

    ## ** extract information
    ## 0-order: param
    param <- coef2(object, as.lava = TRUE)
    n.param <- length(param)
    name.param <- names(param)

    ## 1-order: score
    if(robust){
        score <- score2(object, ssc = ssc)
    }
    
    ## 2-order: variance covariance
    vcov.param <- vcov2(object, ssc = ssc)
    warn <- attr(vcov.param, "warning")
    attr(vcov.param, "warning") <- NULL

    if(robust){
        rvcov.param <- crossprod(iid2(object, ssc = ssc, cluster = cluster))
        hessian <- object$sCorrect$hessian
    }

    ## 3-order: derivative of the variance covariance
    if(!is.null(df)){
        dVcov.param <- object$sCorrect$dVcov.param
        keep.param <- dimnames(dVcov.param)[[3]]
    }
    
    ## ** Computation of the df for the robust case
    if(robust && identical(df, "Satterthwaite") && lava.options()$df.robust==1){

        ## update the score/hessian/derivative at the cluster level
        if(!is.null(cluster)){            
            scoreSave <- score
            hessianSave <- hessian

            score <- matrix(NA, nrow = n.cluster, ncol = NCOL(score),
                            dimnames = list(NULL, colnames(score)))
            hessian <- array(NA, dim = c(NCOL(score), NCOL(score), n.cluster),
                             dimnames = list(colnames(score), colnames(score), NULL))            
            for(iCluster in 1:n.cluster){ ## iCluster <- 1
                score[iCluster,] <- colSums(scoreSave[ls.indexCluster[[iCluster]],,drop=FALSE])
                hessian[,,iCluster] <- apply(hessianSave[,,ls.indexCluster[[iCluster]],drop=FALSE],1:2,sum)
            }
            ## compute derivative
            name.3deriv <- dimnames(dVcov.param)[[3]]
            dRvcov.param <- array(NA, dim = c(n.param,n.param,n.param), dimnames = list(name.param,name.param,name.param))
            for(iP in 1:n.param){ ## iP <- 1
                ## if(name.param[iP] %in% name.3deriv){
                    ## term1 <- dVcov.param[,,name.param[iP]] %*% crossprod(score) %*% vcov.param
                ## }else{
                    ## term1 <- matrix(0, nrow = n.param, ncol = n.param)
                ## }
                ## term2 <- vcov.param %*% hessian[iP,,] %*% score %*% vcov.param
                ## dRvcov.param[,,iP] <- term1 + t(term1) + term2 + t(term2)

                term2 <- vcov.param %*% hessian[iP,,] %*% score %*% vcov.param
                dRvcov.param[,,iP] <- term2 + t(term2)
            }
        }else{
            dRvcov.param <- object$sCorrect$dRvcov.param
        }
    }
    
    ## ** normalize linear hypotheses
    if(!is.matrix(linfct)){
        
        res.C <- .createContrast(linfct, name.param = name.param, add.rowname = TRUE)
        linfct <- res.C$contrast
        if(is.null(rhs)){
            rhs <- res.C$null
        }else{
            if(length(rhs)!=length(res.C$null)){
                stop("Incorrect argument \'rhs\' \n",
                     "Must have length ",length(res.C$null),"\n")
            }
            rhs <- setNames(rhs, names(res.C$null))
        }
        
    }else{
        
        if(is.null(colnames(linfct))){
            stop("Argument \'linfct\' must have column names \n")
        }
        if(NCOL(linfct) != n.param){
            stop("Argument \'linfct\' should be a matrix with ",n.param," columns \n")
        }
        if(any(colnames(linfct) %in% name.param == FALSE)){
            txt <- setdiff(colnames(linfct), name.param)
            stop("Argument \'linfct\' has incorrect column names \n",
                 "invalid name(s): \"",paste(txt, collapse = "\" \""),"\"\n")
        }
        if(any(name.param %in% colnames(linfct) == FALSE)){
            txt <- setdiff(name.param, colnames(linfct))
            stop("Argument \'linfct\' has incorrect column names \n",
                 "missing name(s): \"",paste(txt, collapse = "\" \""),"\"\n")
        }
        ## reorder columns according to coefficients
        linfct <- linfct[,name.param,drop=FALSE]
        if(F.test && any(abs(svd(linfct)$d)<1e-10)){
            stop("Argument \'linfct\' is singular \n")
        }
        if(is.null(rhs)){
            rhs <- setNames(rep(0,NROW(linfct)),rownames(linfct))
        }else if(length(rhs)!=NROW(linfct)){
            stop("The length of argument \'rhs\' must match the number of rows of argument \'linfct' \n")
        }
        if(is.null(rownames(linfct))){
            rownames(linfct) <- .contrast2name(linfct, null = rhs)
            rhs <- setNames(rhs, rownames(linfct))
        }
    }

### ** prepare export
    name.hypo <- rownames(linfct)
    n.hypo <- NROW(linfct)

    df.table <- as.data.frame(matrix(NA, nrow = n.hypo, ncol = 5,
                                     dimnames = list(name.hypo,
                                                     c("estimate","std","statistic","df","p-value"))
                                     ))

    ## ** Univariate Wald test
    C.p <- (linfct %*% param) - rhs
    if(robust){
        C.vcov.C <- linfct %*% rvcov.param %*% t(linfct)
    }else{
        C.vcov.C <- linfct %*% vcov.param %*% t(linfct)
    }
    sd.C.p <- sqrt(diag(C.vcov.C))
    stat.Wald <- C.p/sd.C.p

    ## store
    df.table$estimate <- as.numeric(C.p)
    df.table$std <- as.numeric(sd.C.p)
    df.table$statistic <- as.numeric(stat.Wald)

    ##  degrees of freedom
    if(identical(df,"Satterthwaite") && !is.null(dVcov.param)){

        ## univariate
        if(robust == FALSE){
            df.Wald  <- dfSigma(contrast = linfct,
                                vcov = vcov.param,
                                dVcov = dVcov.param,
                                keep.param = keep.param)
        }else if(robust == TRUE){

                df.Wald <- dfSigma(contrast = linfct,
                                   vcov = vcov.param,
                                   dVcov = dVcov.param,
                                   keep.param = keep.param)
                ## df.Wald  <- dfSigma(contrast = linfct,
                ##                     vcov = rvcov.param,
                ##                     dVcov = dRvcov.param * factor.dRvcov,
                ##                     keep.param = name.param)
                ## df.Wald <- dfSigmaRobust(contrast = linfct,
                ##                          vcov = vcov.param,
                ##                          rvcov = rvcov.param,
                ##                          score = score)
            
        }
    }else{
        df.Wald <- rep(Inf, n.hypo)
        df.F <- Inf
    }

    ## store
    df.table$df <- as.numeric(df.Wald)
    df.table$`p-value` <- as.numeric(2*(1-stats::pt(abs(df.table$statistic), df = df.table$df)))
    
    ## ** Multivariate Wald test
    df.table <- rbind(df.table, global = rep(NA,5))
    error <- NULL
    if(F.test){
        ## statistic
        iC.vcov.C <- try(solve(C.vcov.C), silent = TRUE)
        if(!inherits(iC.vcov.C,"try-error")){
            stat.F <- t(C.p) %*% iC.vcov.C %*% (C.p) / n.hypo

            ## df (independent t statistics)
            if(identical(df,"Satterthwaite")){
                svd.tempo <- eigen(iC.vcov.C)
                D.svd <- diag(svd.tempo$values, nrow = n.hypo, ncol = n.hypo)
                P.svd <- svd.tempo$vectors
     
                C.anova <- sqrt(D.svd) %*% t(P.svd) %*% linfct

                nu_m <- dfSigma(contrast = C.anova,
                                vcov = vcov.param,
                                dVcov = dVcov.param,
                                keep.param = keep.param)
                ## nu_m <- dfSigma(contrast = C.anova,
                ##                 vcov = rvcov.param,
                ##                 dVcov = dRvcov.param * factor.dRvcov,
                ##                 keep.param = keep.param)
                ## nu_m <- dfSigmaRobust(contrast = C.anova,
                ##                       vcov = vcov.param,
                ##                       rvcov = rvcov.param,
                ##                       score = score)
        
                EQ <- sum(nu_m/(nu_m-2))
                df.F <- 2*EQ / (EQ - n.hypo)
            }else{
                df.F <- Inf
            }
            ## store
            df.table["global", "statistic"] <- as.numeric(stat.F)
            df.table["global", "df"] <- df.F
            df.table["global", "p-value"] <- 1 - stats::pf(df.table["global", "statistic"],
                                                           df1 = n.hypo,
                                                           df2 = df.table["global", "df"])
        }else{
            warning("Unable to invert the variance-covariance matrix after application of the contrasts \n")
            error <- iC.vcov.C
        }
    }

    ## ** export
    if(as.lava == TRUE){
        level.inf <- (1-conf.level)/2
        level.sup <- 1-level.inf

        level.inf.label <- paste0(100*level.inf,"%")
        level.sup.label <- paste0(100*level.sup,"%")

        df.estimate <- matrix(NA, nrow = n.hypo, ncol = 5,
                              dimnames = list(name.hypo,c("Estimate", "Std.Err", "df", level.inf.label, level.sup.label)))
        df.estimate[,"Estimate"] <- df.table[name.hypo,"estimate"]
        df.estimate[,"Std.Err"] <- df.table[name.hypo,"std"]
        df.estimate[,"df"] <- df.table[name.hypo,"df"]
        df.estimate[,level.inf.label] <- df.table[name.hypo,"estimate"] + stats::qt(level.inf, df = df.table[name.hypo,"df"]) * df.table[name.hypo,"std"]
        df.estimate[,level.sup.label] <- df.table[name.hypo,"estimate"] + stats::qt(level.sup, df = df.table[name.hypo,"df"]) * df.table[name.hypo,"std"]

        out <- list(statistic = setNames(df.table["global","statistic"],"F-statistic"),
                    parameter = setNames(round(df.table["global","df"],2), paste0("df1 = ",n.hypo,", df2")), ## NOTE: cannot not be change to coefficients because of lava
                    p.value = df.table["global","p-value"],
                    method = c("- Wald test -", "", "Null Hypothesis:", name.hypo),
                    estimate = df.estimate,
                    vcov = C.vcov.C,
                    coef = C.p[,1],
                    null = rhs,
                    cnames = name.hypo                    
                    )
        if(robust){
            colnames(out$estimate)[2] <- "robust SE"
        }        
        attr(out, "B") <- linfct
        class(out) <- "htest"
    }else{
        out <- df.table
        attr(out, "warning") <- warn
        attr(out, "contrast") <- linfct
    }
    attr(out,"error") <- error
    return(out)
}

## * compare2.gls
#' @rdname compare2
#' @export
compare2.gls <- compare2.lm

## * compare2.lme
#' @rdname compare2
#' @export
compare2.lme <- compare2.lm

## * compare2.lvmfit
#' @rdname compare2
#' @export
compare2.lvmfit <- compare2.lm

## * compare2.sCorrect
#' @rdname compare2
compare2.sCorrect <- function(object, linfct = NULL, rhs = NULL,
                              robust = FALSE, cluster = NULL, 
                              ssc = object$sCorrect$ssc$type, df = object$sCorrect$df,
                              as.lava = TRUE, F.test = TRUE, conf.level = 0.95){

    class(object) <- setdiff(class(object),"sCorrect")
    return(compare2(object, linfct = linfct, rhs = rhs,
                    robust = robust, cluster = cluster,
                    ssc = ssc, df = df,
                    as.lava = as.lava, F.test = F.test, conf.level = conf.level))

}

## * compare.sCorrect
#' @rdname compare2
compare.sCorrect <- compare2.sCorrect
    
## * dfSigma
##' @title Degree of Freedom for the Chi-Square Test
##' @description Computation of the degrees of freedom of the chi-squared distribution
##' relative to the model-based variance
##'
##' @param contrast [numeric vector] the linear combination of parameters to test
##' @param vcov [numeric matrix] the variance-covariance matrix of the parameters.
##' @param dVcov [numeric array] the first derivative of the variance-covariance matrix of the parameters.
##' @param keep.param [character vector] the name of the parameters with non-zero first derivative of their variance parameter.
##' 
dfSigma <- function(contrast, vcov, dVcov, keep.param){
    ## iLink <- "LogCau~eta"
    C.vcov.C <- rowSums(contrast %*% vcov * contrast) ## variance matrix of the linear combination
    ## C.vcov.C - vcov[iLink,iLink]

    C.dVcov.C <- sapply(keep.param, function(x){
        rowSums(contrast %*% dVcov[,,x] * contrast)
    })
    ## C.dVcov.C - dVcov[iLink,iLink,]
    numerator <- 2 *(C.vcov.C)^2
    ## numerator - 2*vcov[iLink,iLink]^2
    denom <- rowSums(C.dVcov.C %*% vcov[keep.param,keep.param,drop=FALSE] * C.dVcov.C)
    ## denom - t(dVcov[iLink,iLink,]) %*% vcov[keep.param,keep.param,drop=FALSE] %*% dVcov[iLink,iLink,]
    df <- numerator/denom
    return(df)
}

## * dfSigmaRobust
##' @title Degree of Freedom for the Robust Chi-Square Test
##' @description Computation of the degrees of freedom of the chi-squared distribution
##' relative to the robust-based variance
##'
##' @param contrast [numeric vector] the linear combination of parameters to test
##' @param vcov [numeric matrix] the variance-covariance matrix of the parameters.
##' @param rvcov [numeric matrix] the robust variance-covariance matrix of the parameters.
##' @param score [numeric matrix] the individual score for each parameter.
##'
##' @details When contrast is the identity matrix, this function compute the moments of the sandwich estimator
##' and the degrees of freedom of the approximate t-test as described in (Pan, 2002) section 2 and 3.1.
##'
##' @references
##' Wei Pan and Melanie M. Wall, Small-sample adjustments in using the sandwich variance estiamtor in generalized estimating equations. Statistics in medicine (2002) 21:1429-1441.
##' 
dfSigmaRobust <- function(contrast, vcov, rvcov, score){
    
    ## ** prepare
    n <- NROW(score)

    ## apply contrasts
    vcov.S <- vcov %*% t(contrast)

    ## ** compute moments of rvcov
    ## fast
    E.score2 <- crossprod(score)
    iid.score2 <- lapply(1:n, function(iRow){
        (tcrossprod(score[iRow,]) - E.score2/n)^2
    })

    var.rvcov <- t(vcov.S^2) %*% Reduce("+",iid.score2) %*% vcov.S^2

    ## slow    
    E.rvcov <- contrast %*% rvcov %*% t(contrast)
    ## iid.rvcov <- lapply(1:n, function(iRow){
    ##     (t(vcov.S) %*% tcrossprod(score[iRow,]) %*% vcov.S - E.rvcov/n)^2
    ## })
    ## var.rvcov <- Reduce("+",iid.rvcov)

    ## ** export
    df.rvcov <- (2*E.rvcov^2)/var.rvcov
    
    return(setNames(diag(df.rvcov), rownames(contrast)))
}


##----------------------------------------------------------------------
### compare2.R ends here

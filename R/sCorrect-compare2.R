### compare2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 30 2018 (14:33) 
## Version: 
## Last-Updated: feb 11 2020 (17:21) 
##           By: Brice Ozenne
##     Update #: 773
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
#' @param linfct [matrix or vector of character] the linear hypotheses to be tested. Same as the argument \code{par} of \code{\link{createContrast}}.
#' @param rhs [vector] the right hand side of the linear hypotheses to be tested.
#' @param robust [logical] should the robust standard errors be used instead of the model based standard errors?
#' @param cluster [integer vector] the grouping variable relative to which the observations are iid.
#' @param as.lava [logical] should the output be similar to the one return by \code{lava::compare}?
#' @param F.test [logical] should a joint test be performed?
#' @param conf.level [numeric 0-1] the confidence level of the confidence interval.
#' @param transform [function] function to backtransform the estimates and the associated confidence intervals
#' (e.g. \code{exp} if the outcomes have been log-transformed).
#' @param df [logical] should the degree of freedoms of the Wald statistic be computed using the Satterthwaite correction?
#' Otherwise the degree of freedoms are set to \code{Inf}, i.e. a normal distribution is used instead of a Student's t distribution when computing the p-values.
#' @param ssc [logical] should the standard errors of the coefficients be corrected for small sample bias? Argument passed to \code{sCorrect}.
#' @param ... [logical] arguments passed to lower level methods.
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
    function(object, ...) UseMethod("compare2")

## * compare2.lm
#' @rdname compare2
#' @export
compare2.lm <- function(object, ssc = lava.options()$ssc, df = lava.options()$df, ...){

    object.SSC <- sCorrect(object, ssc = ssc, df = df)
    return(compare2(object.SSC, ...))

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
#' @export
compare2.sCorrect <- function(object, linfct = NULL, rhs = NULL,
                              robust = FALSE, cluster = NULL,
                              as.lava = TRUE, F.test = TRUE,
                              conf.level = 0.95, ...){

    df <- object$sCorrect$df
    
    if(is.null(linfct)){ ## necessary for lava::gof to work
        return(lava::compare(object))
    }
    if(!is.logical(robust)){ 
        stop("Argument \'robust\' should be TRUE or FALSE \n")
    }
    if(robust==FALSE && !is.null(cluster)){
        stop("Argument \'cluster\' must be NULL when argument \'robust\' is FALSE \n")
    }

    ## ** extract information
    ## 0-order: param
    param <- coef2(object, as.lava = TRUE)
    n.param <- length(param)
    name.param <- names(param)

    ## 1-order: score
    if(robust){
        score <- score2(object, cluster = cluster)
    }
    
    ## 2-order: variance covariance
    vcov.param <- vcov2(object)
    warn <- attr(vcov.param, "warning")
    attr(vcov.param, "warning") <- NULL
    if(robust){
        rvcov.param <- crossprod(iid2(object, cluster = cluster))
    }

    ## 3-order: derivative of the variance covariance matrices
    if(identical(df, "Satterthwaite")){
        dVcov.param <- object$sCorrect$dVcov.param
        keep.param <- dimnames(dVcov.param)[[3]]

        if(robust && (lava.options()$df.robust != 1)){

            if(!is.null(cluster)){ ## update derivative according to cluster
                hessian <- hessian2(object, cluster = cluster)
                dRvcov.param <- .dRvcov.param(score = score,
                                              hessian = hessian,
                                              vcov.param = vcov.param,
                                              dVcov.param = dVcov.param,
                                              n.param = n.param,
                                              name.param = name.param)
                                              
            }else{
                dRvcov.param <- object$sCorrect$dRvcov.param
            }
        }
    }

    ## ** normalize linear hypotheses
    if(!is.matrix(linfct)){
        res.C <- createContrast(object, linfct = linfct, add.variance = TRUE, rowname.rhs = FALSE)
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

    n.hypo <- NROW(linfct)
    name.hypo <- rownames(linfct)

    ## ** Univariate Wald test
    ## coefficient (used for F.test and lava export)
    C.p <- linfct %*% param
    C.p.rhs <- C.p - rhs

    ## variance (used for F.test and lava export)
    if(robust){
        C.vcov.C <- linfct %*% rvcov.param %*% t(linfct)
    }else{
        C.vcov.C <- linfct %*% vcov.param %*% t(linfct)
    }

    ## df
    if(identical(df, "Satterthwaite")){
        df.Wald  <- dfSigma(contrast = linfct,
                            score = score,
                            vcov = vcov.param,
                            rvcov = rvcov.param,
                            dVcov = dVcov.param,
                            dRvcov = dRvcov.param,
                            keep.param = keep.param,                            
                            type = if(robust){lava.options()$df.robust}else{1})       
    }else{
        df.Wald <- rep(Inf, n.hypo)
    }
    
    ## ** Multivariate Wald test
    error <- NULL
    if(F.test){
        iC.vcov.C <- try(solve(C.vcov.C), silent = TRUE)
        if(!inherits(iC.vcov.C,"try-error")){
            stat.F <- t(C.p.rhs) %*% iC.vcov.C %*% (C.p.rhs) / n.hypo

            ## df (independent t statistics)
            if(identical(df,"Satterthwaite")){
                svd.tempo <- eigen(iC.vcov.C)
                D.svd <- diag(svd.tempo$values, nrow = n.hypo, ncol = n.hypo)
                P.svd <- svd.tempo$vectors
     
                C.anova <- sqrt(D.svd) %*% t(P.svd) %*% linfct

                nu_m  <- dfSigma(contrast = C.anova,
                                 score = score,
                                 vcov = vcov.param,
                                 rvcov = rvcov.param,
                                 dVcov = dVcov.param,
                                 dRvcov = dRvcov.param,
                                 keep.param = keep.param,                            
                                 type = if(robust){lava.options()$df.robust}else{1})
                EQ <- sum(nu_m/(nu_m-2))
                df.F <- 2*EQ / (EQ - n.hypo)
            }else{
                df.F <- Inf
            }
            ## store
            F.res <- c("statistic" = as.numeric(stat.F),
                       "df" = df.F,
                       "p.value" = 1 - stats::pf(stat.F,
                                                 df1 = n.hypo,
                                                 df2 = df.F)
                       )
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
        df.estimate[,"Estimate"] <- C.p
        df.estimate[,"Std.Err"] <- sqrt(diag(C.vcov.C))
        df.estimate[,"df"] <- df.Wald
        df.estimate[,level.inf.label] <- df.estimate[,"Estimate"] + stats::qt(level.inf, df = df.estimate[,"df"]) * df.estimate[,"Std.Err"]
        df.estimate[,level.sup.label] <- df.estimate[,"Estimate"] + stats::qt(level.sup, df = df.estimate[,"df"]) * df.estimate[,"Std.Err"]

        out <- list(statistic = setNames(F.res["statistic"],"F-statistic"),
                    parameter = setNames(round(F.res["df"],2), paste0("df1 = ",n.hypo,", df2")), ## NOTE: cannot not be change to coefficients because of lava
                    p.value = F.res["p.value"],
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
        if(length(unique(df.Wald))==1){
            df.Wald <- df.Wald[1]
        }
        out <- list(model = object,
                    linfct = linfct,
                    rhs = unname(rhs),
                    coef = param,
                    vcov = if(robust){rvcov.parm}else{vcov.param},
                    df = df.Wald,
                    alternative = "two.sided",
                    type = NULL,
                    robust = robust,
                    ssc = object$sCorrect$ssc$type,
                    global = if(F.test){F.res}else{NULL})
        class(out) <- c("glht2","glht")
    }
    attr(out,"warning") <- warn
    attr(out,"error") <- error
    return(out)

}

## * dfSigma
##' @title Degree of Freedom for the Chi-Square Test
##' @description Computation of the degrees of freedom of the chi-squared distribution
##' relative to the model-based variance
##'
##' @param contrast [numeric vector] the linear combination of parameters to test
##' @param score [numeric matrix] the individual score for each parameter.
##' @param vcov [numeric matrix] the model-based variance-covariance matrix of the parameters.
##' @param rvcov [numeric matrix] the robust variance-covariance matrix of the parameters.
##' @param dVcov [numeric array] the first derivative of the model-based variance-covariance matrix of the parameters.
##' @param dRvcov [numeric array] the first derivative of the robust variance-covariance matrix of the parameters.
##' @param keep.param [character vector] the name of the parameters with non-zero first derivative of their variance parameter.
##' @param type [integer] 1 corresponds to the Satterthwaite approximation of the the degrees of freedom applied to the model-based variance,
##' 2 to the Satterthwaite approximation of the the degrees of freedom applied to the robust variance,
##' 3 to the approximation described in (Pan, 2002) section 2 and 3.1.
##'
##' @references
##' Wei Pan and Melanie M. Wall, Small-sample adjustments in using the sandwich variance estiamtor in generalized estimating equations. Statistics in medicine (2002) 21:1429-1441.
##' 
dfSigma <- function(contrast, score, vcov, rvcov, dVcov, dRvcov, keep.param, type){
    if(type==1){
        C.vcov.C <- rowSums(contrast %*% vcov * contrast) ## variance matrix of the linear combination
        C.dVcov.C <- sapply(keep.param, function(x){
            rowSums(contrast %*% dVcov[,,x] * contrast)
        })
        numerator <- 2 *(C.vcov.C)^2
        denom <- rowSums(C.dVcov.C %*% vcov[keep.param,keep.param,drop=FALSE] * C.dVcov.C)
        df <- numerator/denom
    }else if(type==2){
        C.rvcov.C <- rowSums(contrast %*% rvcov * contrast) ## variance matrix of the linear combination
        C.dRvcov.C <- sapply(keep.param, function(x){
            rowSums(contrast %*% dRvcov[,,x] * contrast)
        })
        numerator <- 2 *(C.rvcov.C)^2
        denom <- rowSums(C.dRvcov.C %*% rvcov[keep.param,keep.param,drop=FALSE] * C.dRvcov.C)
        df <- numerator/denom
    }else if(type==3){
        n <- NROW(score)
        vcov.S <- vcov %*% t(contrast)
        E.score2 <- crossprod(score)
        iid.score2 <- lapply(1:n, function(iRow){
            (tcrossprod(score[iRow,]) - E.score2/n)^2
        })
        var.rvcov <- t(vcov.S^2) %*% Reduce("+",iid.score2) %*% vcov.S^2
        E.rvcov <- contrast %*% rvcov %*% t(contrast)
        df.rvcov <- (2*E.rvcov^2)/var.rvcov
        df <- diag(df.rvcov)
    }
    
    return(setNames(df, rownames(contrast)))
}


##----------------------------------------------------------------------
### compare2.R ends here

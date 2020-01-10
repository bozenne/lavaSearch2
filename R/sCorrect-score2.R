### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: jan 10 2020 (13:26) 
##           By: Brice Ozenne
##     Update #: 2315
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - score2
#' @title  Extract the Individual Score After Small Sample Correction.
#' @description  Extract the individual score from \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} objects.
#' Similar to \code{lava::score} but with small sample correction.
#' @name score2
#'
#' @param object a \code{lm}, \code{gls}, \code{lme}, or \code{lvmfit} object.
#' @param param [optional] the fitted parameters.
#' @param data [optional] the data set.
#' @param bias.correct [logical] should the standard errors of the coefficients be corrected for small sample bias? Only relevant if the \code{sCorrect} function has not yet be applied to the object.
#' @param ... arguments to be passed to \code{sCorrect}.
#'
#' @details If argument \code{p} or \code{data} is not null, then the small sample size correction is recomputed to correct the influence function.
#'
#' @seealso \code{\link{sCorrect}} to obtain \code{lm2}, \code{gls2}, \code{lme2}, or \code{lvmfit2} objects.
#'
#' @return A matrix containing the score relative to each sample (in rows)
#' and each model coefficient (in columns).
#' 
#' @examples
#' n <- 5e1
#' p <- 3
#' X.name <- paste0("X",1:p)
#' link.lvm <- paste0("Y~",X.name)
#' formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))
#'
#' m <- lvm(formula.lvm)
#' distribution(m,~Id) <- sequence.lvm(0)
#' set.seed(10)
#' d <- lava::sim(m,n)
#'
#' ## linear model
#' e.lm <- lm(formula.lvm,data=d)
#' score.tempo <- score2(e.lm, bias.correct = FALSE)
#' colMeans(score.tempo)
#'
#' ## latent variable model
#' e.lvm <- estimate(lvm(formula.lvm),data=d)
#' score.tempo <- score2(e.lvm, bias.correct = FALSE)
#' range(score.tempo-score(e.lvm, indiv = TRUE))
#'
#' @concept small sample inference
#' @export
`score2` <-
  function(object, param, data, ssc, indiv) UseMethod("score2")

## * score2.lm
#' @rdname score2
#' @export
score2.lm <- function(object, param = NULL, data = NULL,
                      ssc = lava.options()$ssc, indiv = FALSE){

    if(is.null(object$sCorrect) || !is.null(param) || !is.null(data) || !identical(object$sCorrect$ssc$type, ssc)){
        object <- sCorrect(object, param = param, data = data, ssc = ssc, df = NA)
    }

    if(indiv){
        return(object$sCorrect$score)
    }else{
        return(colSums(object$sCorrect$score))
    }
}

## * score2.gls
#' @rdname score2
#' @export
score2.gls <- score2.lm

## * score2.lme
#' @rdname score2
#' @export
score2.lme <- score2.lm

## * score2.lvmfit
#' @rdname score2
#' @export
score2.lvmfit <- score2.lm

## * score2.sCorrect
#' @rdname score2
score2.sCorrect <- function(object, param = NULL, data = NULL,
                            ssc = object$sCorrect$ssc$type, indiv = FALSE){
    class(object) <- setdiff(class(object),"sCorrect")
    return(score2(object, param = param, data = data, ssc = ssc, indiv = indiv))
}

## * score.sCorrect
#' @rdname score2
score.sCorrect <- score2.sCorrect

## * .score2
#' @title Compute the Corrected Score.
#' @description Compute the corrected score when there is no missing value.
#' @name score2-internal
#' 
#' @param n.cluster [integer >0] the number of observations.
#' 
#' @keywords internal
.score2 <- function(dmu, dOmega, epsilon, OmegaM1,
                    missing.pattern, unique.pattern, name.pattern,
                    name.param, name.meanparam, name.varparam,
                    n.cluster){
    if(lava.options()$debug){cat(".score2\n")}

    ## ** Prepare
    out.score <- matrix(0, nrow = n.cluster, ncol = length(name.param),
                        dimnames = list(NULL,name.param))
    n.pattern <- length(name.pattern)
    
    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iOmegaM1 <- OmegaM1[[iPattern]]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iEpsilon.OmegaM1 <- epsilon[iIndex,iY,drop=FALSE] %*% iOmegaM1

        ## *** Compute score relative to the mean coefficients
        for(iP in name.meanparam){ # iP <- "Y3~eta"
            out.score[iIndex,iP] <- out.score[iIndex,iP] + rowSums(dmu[[iP]][iIndex,iY,drop=FALSE] * iEpsilon.OmegaM1)
        }
        
        ## *** Compute score relative to the variance-covariance coefficients
        for(iP in name.varparam){ # iP <- "Y2~eta"
            term2 <- - 1/2 * tr(iOmegaM1 %*% dOmega[[iP]][iY,iY,drop=FALSE])            
            term3 <- 1/2 * rowSums(iEpsilon.OmegaM1 %*% dOmega[[iP]][iY,iY,drop=FALSE] * iEpsilon.OmegaM1)
            out.score[iIndex,iP] <- out.score[iIndex,iP] + as.double(term2) + term3
        }        
    }

    ### ** export
    return(out.score)
}


#----------------------------------------------------------------------
### score2.R ends her

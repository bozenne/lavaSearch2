### score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (16:43) 
## Version: 
## last-updated: jan 27 2020 (13:36) 
##           By: Brice Ozenne
##     Update #: 2323
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
#' @param indiv [logical] If \code{TRUE}, the score relative to each observation is returned.
#' Otherwise the total score is returned.
#' @param as.lava [logical] Should the order and the name of the coefficient be the same as those obtained using coef with type = -1.
#' Only relevant for \code{lvmfit} objects.
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
#' eS.lm <- sCorrect(e.lm, ssc = NA, df = NA)
#' score.tempo <- score2(eS.lm)
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
  function(object, indiv, as.lava) UseMethod("score2")

## * score2.sCorrect
#' @rdname score2
score2.sCorrect <- function(object, indiv = FALSE, as.lava = TRUE){
    score <- object$sCorrect$score
    if(as.lava == FALSE){ 
        score <- score[,names(object$sCorrect$skeleton$originalLink2param),drop=FALSE]
        dimames(score) <- list(NULL,as.character(object$sCorrect$skeleton$originalLink2param))
    }

    if(indiv){
        return(score)
    }else{
        return(colSums(score))
    }
}

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

### calcDistMax.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 21 2017 (16:44) 
## Version: 
## last-updated: aug 28 2017 (09:42) 
##           By: Brice Ozenne
##     Update #: 134
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

# {{{ documentation

#' @title Adjust the p.values using the quantiles of the max statistic
#' @description Adjust the p.values using the quantiles of the max statistic.
#' @name calcDistMax
#'
#' @param link the name of the coefficients to test.
#' @param statistic the observed statistic relative to the coefficients to test.
#' @param iid zero-mean iid decomposition of the observed coefficients used to compute the statistic.
#' @param seqSelected the vector of coefficients that have already been selected.
#' @param seqIID zero-mean iid decomposition of the tests to condition on.
#' @param seqQuantile critical threshold of the tests to condition on.
#' If not \code{NULL} the values should correspond the variable in to the first column(s) of the argument iid.
#' @param df the degree of freedom for the t statistic.
#' @param method the method used to compute the p.values. Can be \code{"integration"}, \code{"boot-wild"}, or \code{"boot-norm"}.
#' See the detail section.
#' @param alpha the significance threshold for retaining a new link
#' @param ncpus the number of cpu to use for parellel computations
#' @param initCpus should the cpus be initialized.
#' @param n.sim the total number of simulations.
#' @param n.repMax the maximum number of rejection when using "\code{"boot-wild"} or \code{"boot-norm"}.
#' @param trace should the execution of the function be traced.
#' 
#' @examples 
#' library(mvtnorm)
#'
#' set.seed(10)
#' n <- 100
#' p <- 5
#' link <- letters[1:p]
#' 
#' X.iid <- rmvnorm(n, mean = rep(0,p), sigma = diag(1,p))
#' colnames(X.iid) <- link
#' statistic <- setNames(1:p,link)
#'
#' 
#' \dontrun{
#'  n.sim <- 1e3
#' }
#' \dontshow{
#'  n.sim <- 10
#' }
#' 
#' r1 <- calcDistMax(link, statistic = statistic, iid = X.iid, 
#'             seqIID = NULL, seqQuantile =  NULL, method = "integration",
#'             trace = FALSE, alpha = 0.05, ncpus = 1, initCpus = TRUE, df = 1e6)
#' 
#' r2 <- calcDistMax(link, statistic = statistic, iid = X.iid,
#'             seqIID = NULL, seqQuantile =  NULL, method = "boot-wild",
#'             trace = FALSE, alpha = 0.05, ncpus = 1, initCpus = TRUE, n.sim = n.sim, df = 1e6)
#'
#' r3 <- calcDistMax(link, statistic = statistic, iid = X.iid,
#'             seqIID = NULL, seqQuantile =  NULL, method = "boot-norm",
#'             trace = FALSE, alpha = 0.05, ncpus = 1, initCpus = TRUE, n.sim = n.sim, df = 1e6)
#'
#' r4 <- calcDistMax(link, statistic = statistic, iid = X.iid,
#'             seqIID = NULL, seqQuantile =  NULL, method = "bootstrap",
#'             trace = FALSE, alpha = 0.05, ncpus = 1, initCpus = TRUE, n.sim = n.sim, df = 1e6)
#' 
#' rbind(integration = c(r1$p.adjust, quantile = r1$z),
#'       bootWild    = c(r2$p.adjust, quantile = r2$z),
#'       bootNorm    = c(r3$p.adjust, quantile = r3$z),
#'       boostrap    = c(r4$p.adjust, quantile = r4$z))
#' 
# }}}


#' @rdname calcDistMax
#' @export
calcDistMax <- function(link, statistic, iid, seqIID, seqSelected, seqQuantile, df,
                        method, alpha, ncpus, initCpus, n.sim, trace, n.repMax = 100){

    p.iid <- NCOL(iid)
    n <- NROW(iid)
    conditional <- length(seqIID)>0

    iid.statistic <- scale(iid, center = FALSE, scale = TRUE)        
    Sigma.statistic <- cov(iid.statistic, use = "pairwise.complete.obs")
    
    out <- list(p.adjust = NULL, z = NULL, Sigma = Sigma.statistic)
    
    if(method == "integration"){
        if(trace > 0){ cat("Computation of multivariate normal probabilities to adjust the p.values: ") }
        
        # {{{ warper
        if(conditional==FALSE){
            # {{{ first step
           
            ## compute significance threshold
            out$z <- mvtnorm::qmvt(1-alpha, delta = rep(0,p.iid), sigma = Sigma.statistic, df = df, tail = "both.tails")$quantile
       
            ## adjust p.values
            warperP <- function(name){
                value <- abs(statistic[name])
                if(!is.na(value)){
                    p <- mvtnorm::pmvt(lower = -value, upper = value,
                                       delta = rep(0, p.iid), sigma = Sigma.statistic, df = df)
                    return(1-p)
                }else{
                    return(NA)
                }   
            }
            # }}}          
        }else{

            stop("conditional not implemented \n")
            # {{{ following steps
            n.conditional <- NCOL(seqQuantile)
            seqQuantile.selected <- sapply(1:n.conditional, function(iC){
                seqQuantile[seqSelected[iC],iC]
            })
            
            ## truncation
            vec.lower <- rep(-Inf,p.iid)
            vec.upper <- rep(+Inf,p.iid)
            vecSelected.lower <- rep(-Inf,n.conditional)
            vecSelected.upper <- rep(+Inf,n.conditional)
            vecNselected.lower <- -seqQuantile.selected
            vecNselected.upper <- seqQuantile.selected
            vecAll.lower <- c(vec.lower,vecSelected.lower,vecNselected.lower)
            vecAll.upper <- c(vec.upper,vecSelected.upper,vecNselected.upper)

            pAll.iid <- length(vecAll.upper)
            
            ## limit of integration            
            vecSelected.lowerX <- -seqQuantile.selected
            vecSelected.upperX <- seqQuantile.selected
            vecNselected.lowerX <- rep(-Inf,n.conditional)
            vecNselected.upperX <- rep(+Inf,n.conditional)

            ## influence function
            seqIID.statistic <- lapply(seqIID, scale, center = FALSE, scale = TRUE)
            iidSelected.statistic <- sapply(1:n.conditional, function(iC){
                seqIID.statistic[[iC]][,seqSelected[iC],drop=FALSE]
            })
                
            ## compute significance threshold
            wraperQ <- function(v){
                resV <- 1-tmvtnorm::ptmvnorm(lower = vec.lower, upper = vec.upper,
                                             lowerx = rep(-v,p.iid), upperx = rep(+v,p.iid), 
                                             mean = rep(0,p.iid), sigma = Sigma.statistic)
                return(alpha-resV)
            }

            ## adjust p.values
            warperP <- function(name){ # name <- link[1]
                value <- abs(statistic[name])
                
                if(!is.na(value)){
                    vec.lowerX <-rep(-value,p.iid)
                    vec.upperX <-rep(value,p.iid)
                    vecAll.lowerX <- c(vec.lowerX,vecSelected.lowerX,vecNselected.lowerX)
                    vecAll.upperX <- c(vec.upperX,vecSelected.upperX,vecNselected.upperX)
                    vecAll.lowerX2 <- vecAll.lowerX
                    vecAll.lowerX2[1:p.iid] <- -Inf
                    vecAll.upperX2 <- vecAll.upperX
                    vecAll.upperX2[1:p.iid] <- +Inf
                    
                    iidNSelected.statistic <- sapply(1:n.conditional, function(iC){
                        seqIID.statistic[[iC]][,name]
                    })
                    SigmaAll.statistic <- cov(cbind(iid.statistic,iidSelected.statistic,iidNSelected.statistic),
                                              use = "pairwise.complete.obs")

                    ## P[Tnew>znew|Tselected>zold,Tnselected<zold] < P[Tnew>znew|Tselected>zold,Tnselected<Tselected]
                    ## P[Tnew>znew|Tselected>zold,Tnselected<zold]
                    ## = P[Tnew>znew,Tselected>zold|Tnselected<zold]/P[Tselected>zold|Tnselected<zold]
                    ## = (1-P[Tnew<znew,Tselected<zold|Tnselected<zold])/(1-P[Tselected<zold|Tnselected<zold])
                    p0 <- mvtnorm::pmvt(lower = vec.lowerX, upper = vec.upperX,
                                        delta = rep(0,p.iid), sigma = Sigma.statistic, df = df)
                    
                    p1 <- tmvtnorm::ptmvt(lower = vecAll.lower, upper = vecAll.upper,
                                          lowerx = vecAll.lowerX, upperx = vecAll.upperX,
                                          mean = rep(0,pAll.iid), sigma = SigmaAll.statistic, df = df)
                    p2 <- tmvtnorm::ptmvt(lower = vecAll.lower, upper = vecAll.upper,
                                          lowerx = vecAll.lowerX2, upperx = vecAll.upperX2,
                                          mean = rep(0,pAll.iid), sigma = SigmaAll.statistic)


                    p1 <- mvtnorm::pmvt(lower = c(-0.01,-2), upper = c(0.01,2),
                                        delta = c(0,0), sigma = diag(1,2), df = df)
                    p2 <- mvtnorm::pmvt(lower = c(-Inf,-2), upper = c(Inf,2),
                                        delta = c(0,0), sigma = diag(1,2), df = df)
                    return( (1-p1)/(1-p2) )
                }else{
                    return(NA)
                }   
            }       
            # }}}
        }
        # }}}

        # {{{ parallel computations
        if(initCpus){
            cl <- parallel::makeCluster(ncpus)
            doSNOW::registerDoSNOW(cl)
        }

        value <- NULL # for CRAN check
        out$p.adjust <- foreach::`%dopar%`(
                                     foreach::foreach(value = link,
                                                      .packages = c("tmvtnorm","mvtnorm"),
                                                      .export = "calcDistMax",
                                                      .combine = "c"),
                                     {
                                         warperP(value)
                                     })

        if(initCpus){
            parallel::stopCluster(cl)
        }
    
        if(trace > 0){ cat("done \n") }
        # }}}
        
  }else{

      # {{{ bootstrap
      if(method %in% "boot-norm"){
          Sigma.iid <- cov(iid)
      }else{
          Sigma.iid <- NULL
      }
    
      if(ncpus>1){
          n.simCpus <- rep(round(n.sim/ncpus),ncpus)
          n.simCpus[1] <- n.sim-sum(n.simCpus[-1])
      }else{
          n.simCpus <- n.sim
      }
    
      warper <- function(iid, sigma, n, method){
          if(method == "boot-wild"){
              e <- rnorm(n,mean=0,sd=1)
              iid.sim <- sapply(1:p.iid,function(x){e*iid[,x]})        
          }else if(method == "boot-norm"){
              iid.sim <- MASS::mvrnorm(n,rep(0,p.iid),sigma)            
          }else if(method == "bootstrap"){
              iid.sim <- iid[sample.int(n, replace = TRUE),]
          }
          # apply(iid.sim,2,mean)
          # apply(iid.sim,2,sd)
          Test <- apply(iid.sim,2,function(x){sqrt(n)*mean(x)/sd(x)})                
          return(max(abs(Test)))
      }
    
      if(trace > 0){ cat("simulation to get the 95% quantile of the max statistic: ") }
      
      if(initCpus){
          cl <- parallel::makeCluster(ncpus)
          doSNOW::registerDoSNOW(cl)
      }
      
      i <- NULL # for CRAN check
      distMax <- foreach::`%dopar%`(
                              foreach::foreach(i = 1:ncpus, .packages =  c("MASS"),
                                               .export = "calcDistMax",
                                               .combine = "c"),
                              {
                                  replicate(n.simCpus[i], warper(iid = iid,
                                                                 sigma = Sigma.iid,
                                                                 n = n, method = method))
                              })

      if(initCpus){
          parallel::stopCluster(cl)
      }

      if(trace > 0){ cat("done \n") }
      out$z <- quantile(distMax, probs = 1-alpha)
      out$p.adjust <- sapply(abs(statistic), function(x){mean(distMax>x)})
      # }}}
  }

    return(out)
}


#----------------------------------------------------------------------
### calcDistMax.R ends here

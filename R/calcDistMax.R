### calcDistMax.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: jun 21 2017 (16:44) 
## Version: 
## last-updated: jun 22 2017 (16:30) 
##           By: Brice Ozenne
##     Update #: 10
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
#' @param statistic the observed statistic
#' @param iid zero-mean iid decomposition of the observed coefficients used to compute the statistic.
#' @param mu estimated value for the coefficients
#' @param conditional values to condition on. 
#' If not \code{NULL} the values should correspond the variable in to the first column(s) of the argument iid.
#' @param method the method used to compute the p.values. Can be \code{"integration"}, \code{"boot-wild"}, or \code{"boot-norm"}.
#' See the detail section.
#' @param alpha the significance threshold for retaining a new link
#' @param ncpus the number of cpu to use for parellel computations
#' @param n.sim the total number of simulations.
#' @param n.repMax the maximum number of rejection when using "\code{"boot-wild"} or \code{"boot-norm"}.
#' @param trace should the execution of the function be traced.
#' 
#' @examples 
#' library(mvtnorm)
#' 
#' n <- 100
#' p <- 5
#' X.iid <- rmvnorm(n, mean = rep(0,p), sigma = diag(1,p))
#'
#' calcDistMax(1:p, iid = X.iid, mu = rep(1,p), conditional =  NULL, method = "integration",
#'             trace = FALSE, alpha = 0.05, ncpus = 1)
#' 
#' calcDistMax(1:p, iid = X.iid, mu = rep(1,p), conditional =  NULL, method = "boot-wild",
#'             trace = FALSE, alpha = 0.05, ncpus = 1, n.sim = 10)
#' 
# }}}


#' @rdname calcDistMax
#' @export
calcDistMax <- function(statistic, iid, mu, conditional, method, alpha, ncpus, n.sim, trace, n.repMax = 100){
  
  method <- match.arg(method, c("integration","boot-wild","boot-norm"))
  p.iid <- NCOL(iid)
  n <- NROW(iid)
  
  out <- list(p.adjust = NULL, z = NULL)


  if(method == "integration"){
    if(trace > 0){ cat("Computation of multivariate normal probabilities to adjust the p.values: ") }
      
    iid.statistic <- sapply(1:p.iid, function(col){(mu[col]+iid[,col])/sd(iid[,col])})
    Sigma.statistic <- cov(iid.statistic)
    # mean(iid.statistic)-statistic

    vec.lower <- rep(-Inf,p.iid)
    vec.upper <- rep(+Inf,p.iid)
    
    if(all(conditional==0)){
        ## compute significance threshold
        out$z <- mvtnorm::qmvnorm(p = 1-alpha/2, mean = rep(0,p.iid), sigma = Sigma.statistic)$quantile

        ## adjust p.values
        warperP <- function(value){
            if(!is.na(value)){
                p <- mvtnorm::pmvnorm(lower = -value, 
                                      upper = value, mean = rep(0, p.iid), sigma = Sigma.statistic)
                return(1-p)
            }else{
                return(NA)
            }   
        }
    }else{
        ## compute significance threshold
        wraperQ <- function(v){
            resV <- 1-tmvtnorm::ptmvnorm(lower = vec.lower, upper = vec.upper,
                                         lowerx = rep(-v,p.iid), upperx = rep(+v,p.iid), 
                                         mean = rep(0,p.iid), sigma = Sigma.statistic)
            return(alpha-resV)
        }
        
        max.stat <- 2+log(p.iid)
        out$z <- uniroot(wraperQ, interval = c(1,max.stat))$root
        
        ## adjust p.values
        index.Condn0 <- which(conditional!=0)
        vec.lower[index.Condn0] <- -conditional[index.Condn0]      
        vec.upper[index.Condn0] <- conditional[index.Condn0]
      
        warperP <- function(value){
            if(!is.na(value)){
                p <- tmvtnorm::ptmvnorm(lower = vec.lower, upper = vec.upper,
                                        lowerx = rep(-value,p.iid), upperx = rep(value,p.iid), 
                                        mean = rep(0,p.iid), sigma = Sigma.statistic)
                return(1-p)
            }else{
                return(NA)
            }   
        }       
        
    }

    cl <- parallel::makeCluster(ncpus)
    doSNOW::registerDoSNOW(cl)

    value <- NULL # for CRAN check
    out$p.adjust <- foreach::`%dopar%`(
                                 foreach::foreach(value = abs(statistic),
                                                  .packages = c("tmvtnorm","mvtnorm"),
                                                  .export = "calcDistMax",
                                                  .combine = "c"),
                                 {
                                     warperP(value)
                                 })
    
    parallel::stopCluster(cl)
    
    if(trace > 0){ cat("done \n") }
    
  }else{
    
    if(method %in% "boot-norm"){
      Sigma <- cov(iid)
    }else{
      Sigma <- NULL
    }
    
    if(ncpus>1){
      n.simCpus <- rep(round(n.sim/ncpus),ncpus)
      n.simCpus[1] <- n.sim-sum(n.simCpus[-1])
    }else{
      n.simCpus <- ncpus
    }
    
    warper <- function(iid, mu, sigma, n, method, conditional){
      iRep <- 1
      test.condition <- FALSE
      #    if(!is.null(conditional)){browser()}
      while((iRep <= n.repMax) && (test.condition==FALSE)){
        
        if(method == "boot-wild"){
          e <- rnorm(n,mean=0,sd=1)
          iid.sim <- sapply(1:p.iid,function(x){e*iid[,x]+mu[x]})        
        }else if(method == "boot-norm"){
          iid.sim <- MASS::mvrnorm(n,mu,Sigma)            
        }else {
          stop("method must be \"wild\" or \"normal\"")
        }
        Test <- apply(iid.sim,2,function(x){sqrt(n)*mean(x)/sd(x)})
        if(!is.null(conditional)){
          cmatch <- match(names(conditional),colnames(iid))
          test.condition <- all(abs(Test[cmatch])>conditional)
          Test <- Test[-cmatch]
          iRep <- iRep+1
        }else{
          test.condition <- TRUE
        }
        
      }
      
      return(max(abs(Test)))
    }
    
    if(trace > 0){ cat("simulation to get the 95% quantile of the max statistic: ") }
    conditional.n0 <- conditional[conditional!=0]
    if(length(conditional.n0)==0){
      conditional.n0 <- NULL
    }
    mu[conditional==0] <- 0 # test null hypothesis
    
    cl <- parallel::makeCluster(ncpus)
    doSNOW::registerDoSNOW(cl)

    i <- NULL # for CRAN check
    distMax <- foreach::`%dopar%`(
      foreach::foreach(i = 1:ncpus, .packages =  c("MASS"),
                       .export = "calcDistMax",
                       .combine = "c"),
      {
        replicate(n.simCpus[i], warper(iid = iid,
                                   mu = mu, sigma = Sigma,
                                   n = n, 
                                   method = method, conditional = conditional.n0))
      })
    
    parallel::stopCluster(cl)
    if(trace > 0){ cat("done \n") }
    out$z <- quantile(distMax, probs = 1-alpha)
    out$p.adjust <- sapply(abs(statistic), function(x){mean(distMax>x)})
  }
  
  return(out)
}


#----------------------------------------------------------------------
### calcDistMax.R ends here

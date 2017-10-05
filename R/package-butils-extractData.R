#' @title Extract data from a model
#' 
#' @description Extract data from a model using \code{nlme::getData}, \code{riskRegression::coxDesign} or \code{model.frame}.. 
#' If it fails it will try to extract it by its name according to \code{model$call$data}.
#' 
#' @param object the fitted model.
#' @param model.frame should the data be extracted after transformation (e.g. using model frame)
#' or should the original dataset be extracted.
#' @param convert2dt should the object containing the data be converted into a data.table?
#'  
#' @examples
#' set.seed(10)
#' n <- 101
#'
#' #### linear regression ####
#' Y1 <- rnorm(n, mean = 0)
#' Y2 <- rnorm(n, mean = 0.3)
#' Id <- findInterval(runif(n), seq(0.1,1,0.1))
#' df <- rbind(data.frame(Y=Y1,G="1",Id = Id),
#'            data.frame(Y=Y2,G="2",Id = Id)
#'            )
#' m.lm <- lm(Y ~ G, data = df)
#' extractData(m.lm, model.frame = TRUE)
#' extractData(m.lm, model.frame = FALSE)
#' 
#' library(nlme)
#' m.gls <- gls(Y ~ G, weights = varIdent(form = ~ 1|Id), data = df)
#' extractData(m.gls)
#' m.lme <- lme(Y ~ G, random = ~ 1|Id, data = df)
#' extractData(m.lme)
#' 
#' library(lava)
#' e <- estimate(lvm(Y ~ G), data = df)
#' extractData(e)
#' extractData(e, model.frame = TRUE)
#' 
#' #### survival ####
#' library(riskRegression)
#' library(survival)
#' dt.surv <- sampleData(n, outcome = "survival")
#' m.cox <- coxph(Surv(time, event) ~ X1 + X2, data = dt.surv, x = TRUE, y = TRUE)
#' extractData(m.cox, model.frame = FALSE)
#' extractData(m.cox, model.frame = TRUE)
#' m.cox <- coxph(Surv(time, event) ~ strata(X1) + X2, data = dt.surv, x = TRUE, y = TRUE)
#' extractData(m.cox, model.frame = TRUE)
#' 
#' #### nested fuuctions ####
#' fct1 <- function(m){
#'    fct2(m)
#' }
#' fct2 <- function(m){ 
#'    extractData(m)
#' }
#' fct1(m.gls)
#' @export
extractData <- function(object, model.frame = FALSE, convert2dt = TRUE){

    ## check arguments
    validLogical(convert2dt, valid.length = 1)
    validLogical(model.frame, valid.length = 1)
    
    if(model.frame){ ## use extractors 
        if(any(class(object) %in% c("gls","gnls","lme","lmList","nlme","nls"))){ # nlme package
      
            name.data <- as.character(object$call$data)
      
            # assign the dataset to the object if not in the current environment
            if(name.data %in% ls() == FALSE){
                object$data <- findInParent(name.data, environment())
            }
      
            data <- try(nlme::getData(object), silent = TRUE)
      
        }else if(any(class(object) %in% c("coxph","cph"))){
      
            requireNamespace("riskRegression")
            data <- try(riskRegression::coxDesign(object), silent = TRUE)
            strataVar <- riskRegression::coxVariableName(object)$stratavars.original
      
            if(length(strataVar)>0){ 
        
                data2 <- getInParentEnv(as.character(object$call$data), environment())
        
                data2 <- as.data.table(data2)
                data <- cbind(data, data2[,.SD,.SDcols = strataVar])
        
            }
        }else{
            data <- try(model.frame(object), silent = TRUE)
        }
    
        ## check error
        if("try-error" %in% class(data)){
            stop(data)
        }
    
    }else{
        data <- getInParentEnv(as.character(object$call$data), environment())
    
        if(is.null(data)){
            stop("Could not extract the data from the model \n")
        }  

  }  
    
  ## conversion to data.table
  if(convert2dt){
    if(data.table::is.data.table(data)){
      data <- copy(data)
    }else{
      data <- as.data.table(data)
    }
  }
  
  ## export
  return(data)
}

#' @title Find object in the parent environments
#' 
#' @description Internal function
#' 
#' @param name character string containing the name of the object to get.
#' @param envir the environment from which to look for the object.
findInParent <- function(name, envir){

    frames <- sys.status()
    all.frames <- sapply(1:length(frames$sys.frames), function(x){identical(parent.frame(x),globalenv())})
    index.parents <- which(all.frames==FALSE)
    n.parents <- length(index.parents)
    
    iParent <- 1
    res <- NULL
    while(iParent <= n.parents){ # iParent <- 1
      if(name %in% ls(envir = parent.frame(iParent))){
        res <- get(name, envir = parent.frame(iParent))
        iParent <- n.parents + 1
      }else{
        iParent <- iParent + 1     
      }
    }
    
    return(res)
}

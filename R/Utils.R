### Utils.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 27 2018 (14:32) 
## Version: 
## Last-Updated: nov 25 2019 (15:16) 
##           By: Brice Ozenne
##     Update #: 155
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * formula
## ** selectResponse (Documentation)
#' @title Response Variable of a Formula
#' @description Return the response variable contained in the formula.
#' @name selectResponse
#' 
#' @param object a formula
#' @param format [character] should an object of type call be returned (\code{format = "call"}),
#' or the names of the variables (\code{format = "vars"})
#' @param ... [internal] Only used by the generic method.
#'
#' @return See argument \code{format}.
#' 
#' @examples
#'
#' \dontrun{
#'
#' selectResponse <- lavaSearch2:::selectResponse
#' selectResponse.formula <- lavaSearch2:::selectResponse.formula
#' 
#' selectResponse(Y1~X1+X2)
#' selectResponse(Y1~X1+X2, format = "vars")
#' selectResponse(Surv(event,time)~X1+X2, format = "vars")
#' 
#' selectResponse(Y1~X1+Y1)
#' selectResponse(Y1+Y2~X1+Y1, format = "vars")
#' 
#' selectResponse(~X1+X2)
#' selectResponse(~X1+X2, format = "vars")
#' }
#' 
#' @rdname selectResponse
#' @keywords internal
`selectResponse` <-  function(object, ...) UseMethod("selectResponse")

## ** selectResponse.formula
#' @rdname selectResponse
#' @method selectResponse formula
selectResponse.formula <- function(object, format = "call", ...){
  
  match.arg(format, c("call","vars"))
  
  if(length(object)==3){
    res <- object[[2]]
    if(format == "vars"){
      res <- all.vars(res)
    }
  }else{
    res <- NULL
  }
  
  return(res)
}

## ** selectRegressor (Documentation)
#' @title Regressor of a Formula.
#' @description Return the regressor variables contained in the formula
#' @name selectRegressor
#' 
#' @param object a formula
#' @param format [character] should an object of format call be returned (\code{format = "call"}),
#' or the names of the variables (\code{format = "vars"})
#' @param ... [internal] Only used by the generic method.
#'
#' 
#' @examples
#'
#' \dontrun{
#'
#' selectRegressor <- lavaSearch2:::selectRegressor
#' selectRegressor.formula <- lavaSearch2:::selectRegressor.formula
#' 
#' selectRegressor(Y1~X1+X2)
#' selectRegressor(Y1~X1+X2, format = "vars")
#' 
#' selectRegressor(Y1~X1+Y1)
#' selectRegressor(Y1+Y2~X1+Y1, format = "vars")
#' 
#' selectRegressor(~X1+X2)
#' selectRegressor(~X1+X2, format = "vars")
#' 
#' }
#' @rdname selectRegressor
#' @keywords internal
`selectRegressor` <-  function(object, ...) UseMethod("selectRegressor")

## ** selectRegressor.formula
#' @rdname selectRegressor
#' @method selectRegressor formula
selectRegressor.formula <- function(object, format = "call", ...){
  
  match.arg(format, c("call","vars"))
  
  if(length(object)==3){
    res <- object[[3]]
    
  }else if(length(object)==2){
    res <- object[[2]]
  }else{
    res <- NULL
  }
  if(format == "vars"){
    res <- all.vars(res)
  }
  
  return(res)
}

## ** combineFormula
#' @title Combine formula
#' @description Combine formula by outcome
#' 
#' @param ls.formula a list of formula
#' @param as.formula should a list of formula be returned. Otherwise it will be a list of characters.
#' @param as.unique should regressors appears at most once in the formula
#' 
#' @examples
#' combineFormula(list(Y~X1,Y~X3+X5,Y1~X2))
#' lava.options(symbols = c("~",","))
#' combineFormula(list("Y~X1","Y~X3+X5","Y1~X2"))
#' lava.options(symbols = c("<-","<->"))
#' combineFormula(list("Y<-X1","Y<-X3+X5","Y1<-X2"))
#' 
#' combineFormula(list(Y~X1,Y~X3+X1,Y1~X2))
#' combineFormula(list(Y~X1,Y~X3+X1,Y1~X2), as.formula = FALSE)
#' combineFormula(list(Y~X1,Y~X3+X1,Y1~X2), as.unique = TRUE)
#' 
#' lava.options(symbols = c("~","~~"))
#' combineFormula(list("Y~X1","Y~X3","Y1~X2"))
#' 
#' @export
combineFormula <- function(ls.formula, as.formula = TRUE, as.unique = FALSE){
  
  if(length(ls.formula)==0){return(NULL)}
  ls.Vars <- initVarLinks(ls.formula, format = "list")
  
  ls.endogeneous <- ls.Vars$var1
  ls.X <- ls.Vars$var2
  endogenous <- unique(ls.endogeneous)
  n.endogeneous <- length(endogenous)
  
  ls.formula2 <- vector(n.endogeneous, mode = "list")
  for(iterE in 1:n.endogeneous){
    X <- unlist(ls.X[which(ls.endogeneous==endogenous[iterE])])
    if(as.unique){X <- unique(X)}
    txt <- paste(endogenous[iterE],"~",paste(X, collapse = " + "))
    if(as.formula){ls.formula2[[iterE]] <- as.formula(txt)}else{ls.formula2[[iterE]] <- txt}
  }
  
  return(ls.formula2)
}



## ** formula2character
#' @title formula character conversion
#' @description Conversion of formula into character string or vice versa
#' @name convFormulaCharacter
#' 
#' @param f a formula.
#' @param type should the normal formula operator be used (\code{"formula"}) or the one of lava.option (\code{"symbols"} or \code{"symbol"}).
#' 
#' @examples
#' formula2character(Y1~X1+X2)
#' formula2character(Y1~X1+X2, type = "symbols")

#' @rdname convFormulaCharacter
#' @export
formula2character <- function(f, type = "formula"){
  
  match.arg(type, choices = c("formula", "symbols"))
  
  if(type == "formula"){
    txt <- paste(deparse(f), collapse = "+")
  }else {
    txt <- as.character(f)
    txt[1] <- lava.options()[[type]][1]
    txt <- paste(txt[2],txt[1],txt[3], sep = "")
  }
  
  return(gsub("[[:blank:]]","",txt))
  
}

## * Miscellaneous
## ** .allPermutations
## .allPermutations(1:3)
## .allPermutations(2:3)
.allPermutations <- function(vec){
    X <- lapply(vec, function(x){
        cbind(x, .allPermutations(setdiff(vec, x)))
    })
    return(unname(do.call(rbind,X)))
}
## ** .combination
#' @title Form all Unique Combinations Between two Vectors
#' @description Form all unique combinations between two vectors (removing symmetric combinations).
#' @name combination
#'
#' @param ... [vectors] elements to be combined.
#'
#' @return A matrix, each row being a different combination.
#' 
#' @examples
#' .combination <- lavaSearch2:::.combination
#' 
#' .combination(1,1)
#' .combination(1:2,1:2)
#' .combination(c(1:2,1:2),1:2)
#' 
#' .combination(alpha = 1:2, beta = 3:4)
#'
#' @keywords internal
.combination <- function(...){

    ## ** normalize arguments
    dots <- list(...)
    if(length(dots)!=2){
        stop("can only handle two vectors \n")
    }
    test.null <- unlist(lapply(dots,is.null))    
    if(any(test.null)){
        return(NULL)
    }
    dots <- lapply(dots,unique)

    ## ** form all combinations
    grid <- expand.grid(dots, stringsAsFactors = FALSE) 
    
    ## ** remove combinations (b,a) when (a,b) is already there
    name1 <- paste0(grid[,1],grid[,2])
    name2 <- paste0(grid[,2],grid[,1])

    if(NROW(grid)>1 && any(name1 %in% name2)){ 

        n.grid <- NROW(grid)
        test.duplicated <- c(FALSE,sapply(2:n.grid, function(iG){
            any(name2[iG] %in% name1[1:(iG-1)]) ## find duplicates
        }))

        grid <- grid[test.duplicated==FALSE,]
    }

    ## ** export
    return(grid)        
}


## ** .combinationDF
.combinationDF <- function(data,
                           detail1, detail2,
                           name1, name2){

    detail <- NULL # [:for CRAN check] subset
    
    if(any(detail1 %in% data$detail) && any(detail2 %in% data$detail) ){
        ls.args <- list(subset(data, subset = detail %in% detail1, select = "param", drop = TRUE),
                        subset(data, subset = detail %in% detail2, select = "param", drop = TRUE))
        names(ls.args) <- c(name1,name2)
    
        return(do.call(.combination, args = ls.args))
        
    }else{
        
        return(numeric(0))
        
    }
}

######################################################################
### Utils.R ends here





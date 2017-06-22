# {{{ initVar
#' @title Normalize var1 and var2
#' @name initVar
#' @description Convert var1 and var2 from formula or covariance to character
#' 
#' @param var1 a character indicating the endogeneous variable or a formula
#' @param var2 an optional character indicating the exogeneous variable
#' @param repVar1 should var1 be duplicated to match var2 length. Only active if format = "list".
#' @param format should the name of the variable be return (format = "list"), a vector of character formula ("txt.formula") or a list of formula ("formula")
#' @param Slink the symbol for regression link
#' @param Scov the symbol for covariance link
#' 
#' @details See test file test/testthat/test-Reduce.R for examples
#'
#' @examples
#' initVar_link(y ~ x1)
#' initVar_link("y ~ x1")
#' initVar_link(y ~ x1 + x2)
#' initVar_link("y ~ x1 + x2")
#' initVar_link(y ~ x1 + x2, repVar1 = TRUE)
#' initVar_link(y ~ x1 + x2, repVar1 = TRUE, format = "formula")
#' initVar_link(y ~ x1 + x2, repVar1 = TRUE, format = "txt.formula")
#' initVar_link("y", "x1", format = "formula")
#'
#' initVar_link("y ~ x1:0|1")
#' 
#'
#' initVar_links(y ~ x1)
#' initVar_links("y ~ x1")
#' initVar_links(c("y ~ x1","y~ x2"))
#' initVar_links(c(y ~ x1,y ~ x2))
#' initVar_links(c("y ~ x1","y~ x2"), format = "formula")
#' initVar_links(c(y ~ x1,y ~ x2), format = "formula")
#' initVar_links(c("y ~ x1","y~ x2"), format = "txt.formula")
#' initVar_links(c(y ~ x1,y ~ x2), format = "txt.formula")

#' @rdname initVar
#' @export
initVar_link <- function(var1, var2, repVar1 = FALSE, format = "list",
                         Slink = lava.options()$symbols[1],
                         Scov = lava.options()$symbols[2]){

    format <- match.arg(format, c("list","txt.formula","formula"))
    
    if(missing(var2)){
        if(class(var1) == "formula"){
            var2 <- select.regressor(var1, type = "vars")
            var1 <- select.response(var1, type = "vars")
            sep <- if(format == "formula"){"~"}else{Slink}
        }else if(grepl(Scov,var1,fixed=TRUE)==TRUE){ ## covariance
            varSplit <- strsplit(var1, split = Scov)[[1]]
            var1 <- trimws(varSplit[1])
            var2 <- trimws(varSplit[2])
            sep <- if(format == "formula"){"~"}else{Scov}
        } else if(grepl(Slink,var1,fixed=TRUE)==TRUE){ ## regression
            varSplit <- strsplit(var1, split = Slink)[[1]]
            var1 <- trimws(varSplit[1])
            var2 <- trimws(varSplit[2])
            sep <- if(format == "formula"){"~"}else{Slink}
        } else {
            var1 <- var1
            var2 <- NA
        }
    }else{

        if(!is.character(var1) || !is.character(var2)){
            stop("\'var1\' and \'var2\' must be characters when both are specified \n")
        }
        sep <- if(format == "formula"){"~"}else{Slink}        
    }
  
  
  #### convert to format
  if(format == "formula"){
    n.var2 <- length(var2)
    var1 <- rep(var1, times = n.var2)
    res <- sapply(1:n.var2, function(i){stats::as.formula(paste(var1[i], var2[i], sep = sep))})
    
  }else if(format == "txt.formula"){
    n.var2 <- length(var2)
    var1 <- rep(var1, times = n.var2)
    res <- sapply(1:n.var2, function(i){paste(var1[i], var2[i], sep = sep)})
    
  }else if(format == "list"){
    if(repVar1 && !missing(var1)){var1 <- rep(var1, length(var2))}
    res <- list(var1 = var1,
                var2 = if(!missing(var2)){var2}else{NULL} 
    )
  }
 
  ## export 
  return(res)
}
# }}}

# {{{ initVar_links
#' @rdname initVar
#' @export
initVar_links <- function(var1, format = "list"){
        
    if("formula" %in% class(var1)){
        res <- initVar_link(var1, repVar1 = TRUE, format = format)
    }else {
        res <- sapply(var1, function(x){
            initVar_link(x, repVar1 = TRUE, format = format)
        })
        if(format == "list"){
            res <- list(var1 = unname(unlist(res["var1",])),
                        var2 = unname(unlist(res["var2",])))
        }else{
            res <- unname(unlist(res))
        }
    }

    return(res)
    
}
# }}}

# {{{ combine.formula
#' @title Combine formula
#' @description Combine formula by outcome
#' 
#' @param ls.formula a list of formula
#' @param as.formula should a list of formula be returned. Otherwise it will be a list of characters.
#' @param as.unique should regressors appears at most once in the formula
#' 
#' @examples
#' combine.formula(list(Y~X1,Y~X3+X5,Y1~X2))
#' lava.options(symbols = c("~",","))
#' combine.formula(list("Y~X1","Y~X3+X5","Y1~X2"))
#' lava.options(symbols = c("<-","<->"))
#' combine.formula(list("Y<-X1","Y<-X3+X5","Y1<-X2"))
#' 
#' combine.formula(list(Y~X1,Y~X3+X1,Y1~X2))
#' combine.formula(list(Y~X1,Y~X3+X1,Y1~X2), as.formula = FALSE)
#' combine.formula(list(Y~X1,Y~X3+X1,Y1~X2), as.unique = TRUE)
#' 
#' lava.options(symbols = c("~","~~"))
#' combine.formula(list("Y~X1","Y~X3","Y1~X2"))
#' 
#' @export
combine.formula <- function(ls.formula, as.formula = TRUE, as.unique = FALSE){
  
  if(length(ls.formula)==0){return(NULL)}
  ls.Vars <- initVar_links(ls.formula, format = "list")
  
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
# }}}

# {{{ select.response
#' @title Response variable of a formula
#' @description Return the reponse variable contained in the formula
#' @name select.response
#' 
#' @param x a formula
#' @param type either return an object of type call (\code{"call"}) or the names of the variables (\code{"vars"})
#' @param ... additional arguments to be passed to the low level functions
#' 
#' @examples 
#' select.response(Y1~X1+X2)
#' select.response(Y1~X1+X2, type = "vars")
#' 
#' select.response(Y1~X1+Y1)
#' select.response(Y1+Y2~X1+Y1, type = "vars")
#' 
#' select.response(~X1+X2)
#' select.response(~X1+X2, type = "vars")

#' @rdname select.response
#' @export
`select.response` <-  function(x,...) UseMethod("select.response")

#' @rdname select.response
#' @export
select.response.formula <- function(x, type = "call", ...){
  
  match.arg(type, c("call","vars"))
  
  if(length(x)==3){
    res <- x[[2]]
    if(type == "vars"){
      res <- all.vars(res)
    }
  }else{
    res <- NULL
  }
  
  return(res)
}
# }}}

# {{{ select.regressor
#' @title Regressor of a formula
#' @description Return the regressor variables contained in the formula
#' @name select.regressor
#' 
#' @param x a formula
#' @param type either return an object of type call (\code{"call"}) or the names of the variables (\code{"vars"})
#' @param ... additional arguments to be passed to the low level functions
#'   
#' @examples 
#' select.regressor(Y1~X1+X2)
#' select.regressor(Y1~X1+X2, type = "vars")
#' 
#' select.regressor(Y1~X1+Y1)
#' select.regressor(Y1+Y2~X1+Y1, type = "vars")
#' 
#' select.regressor(~X1+X2)
#' select.regressor(~X1+X2, type = "vars")


#' @rdname select.regressor
#' @export
`select.regressor` <-  function(x,...) UseMethod("select.regressor")

#' @rdname select.regressor
#' @export
select.regressor.formula <- function(x, type = "call", ...){
  
  match.arg(type, c("call","vars"))
  
  if(length(x)==3){
    res <- x[[3]]
    
  }else if(length(x)==2){
    res <- x[[2]]
  }else{
    res <- NULL
  }
  if(type == "vars"){
    res <- all.vars(res)
  }
  
  return(res)
}
# }}}


#### NOT USED ####

# {{{ formula2character
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
# }}}

# {{{ getS3methodParent 
#' @title Find methods belonging to other classes
#' @description Find methods that apply to the object for each class of the object
#' 
#' @param object an object
#' @param FUN the name of a generic method
#' @param class the level of class
#' @param export should all methods be output (\code{"all"}) or only the first one (\code{"first"})
#' @param value should the method be returned. Else its name will be returned.
#' 
#' @examples
#' m <- lvm.reduced(Y~X)
#' getS3methodParent(m, "coef")
#' getS3methodParent(m, "coef", value = TRUE)
#' 
#' @export
getS3methodParent <- function(object, FUN, class = 1, export = "all", value = FALSE){
  
  match.arg(export, c("first","all"))
  
   object.class <- class(object)
   n.class <- length(object.class)
   if(n.class == 0){
     stop("\'object\' inherit from no class \n")
   }
   
   
   if(is.character(class)){
     if(class %in% object.class == FALSE){
       stop("\'object\' does not inherit of  class ",class," \n")
     }else{
       class <- object.class[-(1:which(object.class == class))]
     }
   }else{
     if(class %in% 0:(n.class-1) == FALSE){
       stop("\'object\' inherit from ",n.class," class(es)\n",
            "\'class\' cannot be ",class,"\n")
     }else{
       class <- object.class[-(1:class)]
     }
   }
   ls.method <- sapply(class, function(c){getS3method(f = FUN, class = c, optional = TRUE)})
   test <- unlist(lapply(ls.method, length))
   if(all(test == 0)){
     return(NULL)
   }else if(any(test == 0)){
     ls.method[test == 0] <- NULL
   }
   
   if(value == FALSE){
     available.method <- paste(FUN, class[test>0], sep = ".")
     if(export == "first"){
       return(available.method[1])
     }else if(export == "all"){
       return(available.method)
     }
   }else if(export == "first"){
     return(ls.method[[1]])
   }else if(export == "all"){
     return(ls.method)
   }
}
# }}}

# {{{ callS3methodParent
#' @title Call the first method inhereted by the object
#' @description Call the first method inhereted by the object
#' 
#' @param object an object
#' @param FUN the name of a generic method
#' @param class the level of class
#' @param ... additional arguments to be passed to the method that will be called
#' 
#' @examples
#' m <- lvm.reduced(Y~X)
#' callS3methodParent(m, "coef")
#' callS3methodParent(m, FUN = "coef", class = "lvm.reduced")
#' 
#' @export
callS3methodParent <- function(object, FUN, class = 1,...){
  
  return(getS3methodParent(object, FUN = FUN, class = class, export = "first", value = TRUE)(object, ...))
  
}
# }}}

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


# {{{ select.response
#' @title Response variable of a formula
#' @description Return the reponse variable contained in the formula
#' @name select.response
#' 
#' @param x a formula
#' @param type either return an object of type call (\code{"call"}) or the names of the variables (\code{"vars"})
#' @param ... additional arguments to be passed to lower levels functions.
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
`select.response` <-  function(x, ...) UseMethod("select.response")

#' @rdname select.response
#' @method select.response formula
select.response.formula <- function(x, type = "call"){
  
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
#' @param ... additional arguments to be passed to lower levels functions.
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
`select.regressor` <-  function(x, ...) UseMethod("select.regressor")

#' @rdname select.regressor
#' @method select.regressor formula
select.regressor.formula <- function(x, type = "call"){
  
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



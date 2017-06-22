### Lava_modelsearchLR.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 30 2017 (17:58) 
## Version: 
## last-updated: jun 22 2017 (16:20) 
##           By: Brice Ozenne
##     Update #: 27
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Model searching using a likelihood ratio test
#' @description Model searching using a likelihood ratio test
#' 
#' @name modelsearchLR
#' 
#' @return a list
#'
#' @keywords internal
#'
modelsearchLR <- function (x, data = NULL, restricted, link, directive, 
                           ls.LVMargs = NULL, method.p.adjust = "bonferroni",                           
                           display.warnings = TRUE, trace = 1){

     statistic <- p.value <- adjusted.p.value <- NULL # for CRAN check

    # {{{ initialisation
    n.link <- length(link)
    dt.test <- data.table("link" = link,
                          "statistic" = as.numeric(rep(NA,n.link)),
                          "p.value" = as.numeric(rep(NA,n.link)),
                          "adjusted.p.value" = as.numeric(rep(NA,n.link))
                          )

    if(!is.null(ls.LVMargs)){
        add.args <- setdiff(names(x$call), c("","x","data","control"))
        ls.LVMargs <- lapply(add.args, function(arg){x$call[[arg]]})
        names(ls.LVMargs) <- add.args
        if(is.null(data)){
            ls.LVMargs$data <- x$data$model.frame
        }else{
            ls.LVMargs$data <- data
        }
        ls.LVMargs$control <- x$control
        ls.LVMargs$control$trace <- FALSE
    }

    best.test <- -Inf
    best.model <- NULL
    
    # }}}
    
    if(trace > 0){pb <- utils::txtProgressBar(max = n.link, style = 3) }
  
    for (iterI in 1:n.link) {
    
        ls.LVMargs$x <- addLink(x$model, var1 = restricted[iterI,1], var2 = restricted[iterI,2],
                                covariance = (1-directive[iterI]))
        
        if(display.warnings){
            newfit <- tryCatch(do.call(estimate, args = ls.LVMargs),
                               error = function(x){NA},
                               finally = function(x){x})
        }else{
            suppressWarnings(
                newfit <- tryCatch(do.call(estimate, args = ls.LVMargs),
                                   error = function(x){NA},
                                   finally = function(x){x})
            )
        }
    
        if("lvmfit" %in% class(newfit)){ # test lvmfit is not an error
            if(newfit$opt$convergence == 0){ # test whether lvmfit has correctly converged
                compareT <- lava::compare(x,newfit)
                dt.test[iterI,`statistic` := compareT$statistic[[1]]]
                dt.test[iterI,`p.value` := compareT$p.value[[1]]]

                if(dt.test[iterI,`statistic`]>best.test){
                    best.test <- dt.test[iterI,`statistic`]
                    best.model <- newfit
                }
            }
        }
    
        if(trace > 0){ utils::setTxtProgressBar(pb, value = iterI) }    
    }
    if(trace > 0){  close(pb) }

    dt.test[, adjusted.p.value := p.adjust(p.value, method = method.p.adjust)]
    
    #### export 
    return(list(dt.test = dt.test,
                best.test = best.test,
                best.model = best.model))
}

#----------------------------------------------------------------------
### Lava_modelsearchLR.R ends here

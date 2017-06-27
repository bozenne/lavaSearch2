### Lava_modelsearchLR.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 30 2017 (17:58) 
## Version: 
## last-updated: jun 26 2017 (16:08) 
##           By: Brice Ozenne
##     Update #: 46
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
modelsearchLR <- function (x, data, restricted, link, directive, 
                           update.FCT, update.args,
                           method.p.adjust, display.warnings, trace){

    statistic <- p.value <- adjusted.p.value <- NULL # for CRAN check

    # {{{ initialisation
    n.link <- length(link)
    dt.test <- data.table("link" = link,
                          "statistic" = as.numeric(rep(NA,n.link)),
                          "p.value" = as.numeric(rep(NA,n.link)),
                          "adjusted.p.value" = as.numeric(rep(NA,n.link))
                          )

    best.test <- -Inf
    best.model <- NULL
    
    # }}}
    
    if(trace > 0){pb <- utils::txtProgressBar(max = n.link, style = 3) }

    for (iterI in 1:n.link) { # iterI <- 1
        newfit <- update.FCT(x, args = update.args, restricted = restricted[iterI,], directive = directive[iterI])

        if(class(newfit) != "try-error" && !is.na(logLik(newfit))){ 

            if(newfit$opt$convergence == 0 ){ # test whether the model has correctly converged
                if(class(newfit) == "lvmfit"){
                    compareT <- lava::compare(x,newfit)
                    dt.test[iterI,`statistic` := compareT$statistic[[1]]]
                    dt.test[iterI,`p.value` := compareT$p.value[[1]]]
                }else{
                    compareT <- anova(x, newfit)
                    dt.test[iterI,`statistic` := compareT$F[2]]
                    dt.test[iterI,`p.value` := compareT$`Pr(>F)`[2]]
                }
            }else{
                dt.test[iterI,`statistic` := NA]
                dt.test[iterI,`p.value` := NA]
            }
 
            if(!is.na(dt.test[iterI,`statistic`]) && dt.test[iterI,`statistic`]>best.test){
                best.test <- dt.test[iterI,`statistic`]
                best.model <- newfit
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

### Lava_modelsearchLR.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 30 2017 (17:58) 
## Version: 
## last-updated: sep 25 2017 (11:21) 
##           By: Brice Ozenne
##     Update #: 64
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

    # {{{ initialisation
    n.link <- length(link)
    dt.test <- data.table("link" = link,
                          "statistic" = as.numeric(rep(NA,n.link)),
                          "p.value" = as.numeric(rep(NA,n.link)),
                          "adjusted.p.value" = as.numeric(rep(NA,n.link)),
                          "convergence" = as.numeric(rep(NA,n.link)),
                          "coefBeta" = as.numeric(rep(NA,n.link)),
                          "quantile" = as.numeric(rep(NA,n.link))
                          )

    best.test <- -Inf
    best.model <- NULL
    
    # }}}
    
    if(trace > 0){pb <- utils::txtProgressBar(max = n.link, style = 3) }

    for (iterI in 1:n.link) { # iterI <- 1
        newfit <- update.FCT(x, args = update.args, restricted = restricted[iterI,], directive = directive[iterI])

        if(class(newfit) != "try-error" && !is.na(logLik(newfit))){ 

            if(newfit$opt$convergence == 0 ){ # test whether the model has correctly converged
                newCoef.tempo <- coef(newfit)[setdiff(names(coef(newfit)),names(coef(x)))]
                dt.test[iterI, c("coefBeta") := newCoef.tempo]
                if(class(newfit) == "lvmfit"){
                    compareT <- lava::compare(x,newfit)
                    dt.test[iterI, c("statistic") := compareT$statistic[[1]]]
                    dt.test[iterI, c("p.value") := compareT$p.value[[1]]]
                }else{
                    compareT <- anova(x, newfit)
                    dt.test[iterI, c("statistic") := compareT$F[2]]
                    dt.test[iterI, c("p.value") := compareT$`Pr(>F)`[2]]
                }
                dt.test[iterI, c("convergence") := 0]
            }else{
                dt.test[iterI, c("convergence") := 1]                
            }
 
            if(!is.na(dt.test[iterI][["statistic"]]) && dt.test[iterI][["statistic"]]>best.test){
                best.test <- dt.test[iterI][["statistic"]]
                best.model <- newfit
            }
        }    
    
        if(trace > 0){ utils::setTxtProgressBar(pb, value = iterI) }    
    }
    if(trace > 0){  close(pb) }
    dt.test[, c("adjusted.p.value") := p.adjust(.SD$p.value, method = method.p.adjust)]
    
    #### export 
    return(list(dt.test = dt.test,
                best.test = best.test,
                best.model = best.model))
}

#----------------------------------------------------------------------
### Lava_modelsearchLR.R ends here

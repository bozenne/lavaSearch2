### endogenous.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (13:50) 
## Version: 
## Last-Updated: dec 17 2019 (10:49) 
##           By: Brice Ozenne
##     Update #: 85
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

#' @examples
#' #### simulate data ####
#' set.seed(10)
#' dW <- sampleRepeated(10, format = "wide")
#' set.seed(10)
#' dL <- sampleRepeated(10, format = "long")
#' dL$time2 <- paste0("visit",dL$time)
#' 
#' #### linear model ####
#' e.lm <- lm(Y1~X1, data = dW)
#' endogenous(e.lm)
#'
#' #### gls model ####
#' e.gls1 <- gls(Y1~X1, data = dW)
#' endogenous(e.gls1)
#' 
#' e.gls2 <- gls(Y~X1, correlation = corCompSymm(form=~1|id), data = dL)
#' endogenous(e.gls2)
#'
#' e.gls3 <- gls(Y~X1, correlation = corSymm(form=~time|id), data = dL)
#' endogenous(e.gls3)
#'
#' e.gls4 <- gls(Y~X1, weight = varIdent(form=~1|time2), data = dL)
#' endogenous(e.gls4)
#'
#' e.gls5 <- gls(Y~X1, weight = varIdent(form=~1|time2),
#'               correlation = corSymm(form=~time|id), data = dL)
#' endogenous(e.gls5)
#'
#' #### lme model ####
#' e.lme1 <- lme(Y~X1, random = ~1|id, data = dL)
#' endogenous(e.lme1)
#' 
#' e.lme2 <- lme(Y~X1, random = ~1|id, data = dL,
#'               weight = varIdent(form=~1|time2))
#' endogenous(e.lme2)
#' 
endogenous.lm <- function(x, ...){
    return(selectResponse(formula(x), format = "vars"))
}
endogenous.gls <- function(x, format = "wide", cluster = NULL, ...){
    
    ## ** extract information from object
    name.Y <- selectResponse(formula(x), format = "vars")

    class.re <- class(x$modelStruct$reStruct)
    if("NULL" %in% class.re == FALSE){
        varIndex.re <- unique(names(attr(x$modelStruct$reStruct,"plen")))        
    }else{
        varIndex.re <- NULL
    }

    class.cor <- class(x$modelStruct$corStruct)
    if("NULL" %in% class.cor == FALSE){
        formula.cor <- formula(x$modelStruct$corStruct)
        varIndex.cor <- all.vars(nlme::getCovariateFormula(formula.cor))
    }else{
        varIndex.cor <- NULL
    }

    class.var <- class(x$modelStruct$varStruct)
    if("NULL" %in% class.var == FALSE){
        formula.var <- formula(x$modelStruct$varStruct)
        varIndex.var <- all.vars(nlme::getCovariateFormula(formula.var))
    }else{
        varIndex.var <- NULL
    }

    n <- stats::nobs(x)

    ## ** check arguments
    format <- match.arg(format, choices = c("wide","long"))
    
    if(!is.null(cluster) && length(cluster) != n){
        stop("Argument \'cluster\' has wrong length \n",
             "Should be of length ",n," \n")
    }
    if(length(varIndex.re)>1){
        stop("Can only handle reStruct with one grouping variable for the random effect\n")
    }
    validClass.cor <- c("NULL","corCompSymm","corSymm","corStruct")
    if(any(class.cor %in% validClass.cor == FALSE)){
        stop("Can only handle corStruct of class \"corCompSymm\" or \"corSymm\"\n")
    }

    validClass.var <- c("NULL","varIdent","varFunc")
    if(any(class.var %in% validClass.var == FALSE)){
        stop("Can only handle varStruct of class \"varIdent\"\n")
    }

    if(length(varIndex.cor) > 0 && length(varIndex.cor) > 0 && !identical(varIndex.cor,varIndex.cor)){
        stop("Inconsistency between the left hand side of the formula in corStruct and the left hand side of the formula in varStruct. \n",
             "it should be something like: correlation = corStruct(form = ~index|groupA) \n",
             "                             weight = varStruct(form = ~index|groupB) \n")
    }
    ## ** Identify the index and name of the endogenous variables
    if(format == "long" || ("NULL" %in% class.var && "NULL" %in% class.cor && "NULL" %in% class.re)){ ## basic lme models or lm-ish models, order of the lines in the dataset does not matter
        return(name.Y)
    }else if("NULL" %in% class.var && ("NULL" %in% class.cor == FALSE) && length(varIndex.cor)==0){ ## compound symmetry structure
        name.rep <- sort(unique(unlist(attr(x$modelStruct$corStruct, "covariate"))))
        return(paste0(name.Y, name.rep))
    }else if("NULL" %in% class.var && ("NULL" %in% class.cor == FALSE) && length(varIndex.cor)>0){ ## unstructured
        groupValue.cor <- sort(unique(extractData(x)[[varIndex.cor]]))
        return(paste0(name.Y, groupValue.cor))
    }else if("NULL" %in% class.var){ ## random effect
        n.rep <- max(table(x$groups[[varIndex.re]]))
        ls.relevel <- lapply(x$modelStruct$reStruct, function(iR){
            attr(iR,"Dimnames")[[1]]
        })
        if(length(ls.relevel)==1){
            name.Ywide <- ls.relevel[[1]]
        }else{
            name.Ywide <- levels(interaction(do.call(expand.grid, ls.relevel)))
        }
        if(n.rep==length(name.Ywide)){
            return(paste0(name.Y, name.Ywide))
        }else{
            return(paste0(name.Y, 1:n.rep))
        }
    }else{ ## heterochedastic model
        groupValue.var <- attr(x$modelStruct$varStruct,"groupNames")
        return(paste0(name.Y, groupValue.var))
    }

}
endogenous.lme <- endogenous.gls


######################################################################
### endogenous.R ends here

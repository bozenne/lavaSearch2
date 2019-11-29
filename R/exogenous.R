### exogenous.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 18 2019 (13:52) 
## Version: 
## Last-Updated: nov 25 2019 (09:21) 
##           By: Brice Ozenne
##     Update #: 15
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

exogenous.lm <- function(x,...){
    out <- selectRegressor(formula(x), format = "vars")
    name.coef <- names(stats::coef(x))
    if(length(name.coef)>0){
        name.coef <- name.coef[!grepl(":",name.coef)]
    }
    if(length(name.coef)>0){
        name.coef <- name.coef[name.coef!="(Intercept)"]
    }
    attr(out,"levels") <- name.coef
    return(out)
}
exogenous.gls <- exogenous.lm
exogenous.lme <- exogenous.lm

######################################################################
### exogenous.R ends here

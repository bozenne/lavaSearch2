### test1-sCorrect-adjustedResiduals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  7 2018 (12:21) 
## Version: 
## Last-Updated: jan 17 2022 (14:11) 
##           By: Brice Ozenne
##     Update #: 32
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
## rm(list = ls())
if(FALSE){ ## already called in test-all.R
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))
context("sCorrect (residuals)")

## * simulation
n <- 5e1
mSim <- lvm(c(Y1~eta1,Y2~eta1+X2,Y3~eta1+X1,
              Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(mSim) <- eta1~X1+Gender
latent(mSim) <- ~eta1+eta2
categorical(mSim, labels = c("Male","Female")) <- ~Gender
transform(mSim, Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
d <- lava::sim(mSim, n = n, latent = FALSE)
dL <- reshape2::melt(d, id.vars = c("Id","X1","X2","X3","Gender"),
                     measure.vars = c("Y1","Y2","Y3","Z1","Z2","Z3"))
dLred <- dL[dL$variable %in% c("Y1","Y2","Y3"),]

## * linear regression [lvm]
## ** model fit and sCorrect
e.lvm <- estimate(lvm(Y1~X1+X2+Gender), data = d)
e.lm <- lm(Y1~X1+X2+Gender, data = d)
e.gls <- gls(Y1~X1+X2+Gender, data = d, method = "ML")

e2.lvm <- estimate2(e.lvm)


## ** test adjusted residuals
test_that("residuals2 match residuals.lm (lm adjusted)", {
    epsilon.lm <- residuals(e.lm)
    ## X <- model.matrix(e.lm, d)
    ## iH <- diag(1,n,n) - X %*% solve(t(X) %*% X) %*% t(X)
    ## GS1 <- epsilon.lm/diag(iH)^(1/2)

    ## match expectation
    expect_equal(mean(e2.lvm$sCorrect$residuals), 0)

    ## match variance
    expect_equal(mean(e2.lvm$sCorrect$residuals^2), sigma(e.lm)^2,
                 tol = 1e-8)
    
    ## mean(GS1^2)
    ## mean(residuals(e.lm, type = "pearson")^2)
    ## mean(e2.lm$sCorrect$residuals^2)
    ## sigma(e2.lm)^2

    ## match individual components
    factor <- sqrt(NROW(epsilon.lm)/(NROW(epsilon.lm) - length(coef(e.lm))))
    expect_equal(as.double(e2.lvm$sCorrect$residuals),
                 as.double(epsilon.lm)*factor)
})

## * multivariate linear models
## ** model fit
ls.lm <- list(lm(Y1~X1,d),lm(Y2~X2,d),lm(Y3~X1+X3,d))
e.lvm <- estimate(lvm(Y1~X1,Y2~X2,Y3~X1+X3), data = d)

e2.lvm <- estimate2(e.lvm)

test_that("residuals2 match residuals.lm", {

    epsilon.lm <- lapply(ls.lm, residuals)
    ## X <- lapply(ls.lm, model.matrix)
    ## iH <- lapply(X, function(x){
    ##     diag(1,NROW(x),NROW(x)) - x %*% solve(t(x) %*% x) %*% t(x)
    ## })
    ## GS <- mapply(iH,epsilon.lm, FUN = function(x,y){
    ##     y/diag(x)^(1/2)
    ## })

    ## match expectation
    expect_equal(as.double(colMeans(e2.lvm$sCorrect$residuals)),
                 c(0,0,0),
                 tol = 1e-8)
    
    ## match variance
    expect_equal(unname(colMeans(e2.lvm$sCorrect$residuals^2)),
                 sapply(ls.lm,sigma)^2,
                 tol = 1e-8)

    ## match individual components
    ls.GS <- lapply(ls.lm, function(x){
        iEpsilon <- residuals(x)
        factor <- sqrt(NROW(iEpsilon)/(NROW(iEpsilon) - length(coef(x))))
        return(iEpsilon*factor)
    })

    
    expect_equal(unname(e2.lvm$sCorrect$residuals),
                 unname(do.call(cbind,ls.GS)))
})

## * mixed model: CS [lvm]
m <- lvm(c(Y1[0:sigma]~1*eta,
           Y2[0:sigma]~1*eta,
           Y3[0:sigma]~1*eta,
           eta~X1+X2))
latent(m) <- ~eta
e.lvm <- estimate(m, d)

e.lme <- nlme::lme(value ~ X1 + X2,
                   random =~1| Id,
                   data = dLred, method = "ML")

e.gls <- nlme::gls(value ~ X1 + X2,
                   correlation = corCompSymm(form = ~1| Id),
                   data = dLred, method = "ML")

e2.lvm <- estimate2(e.lvm)

test_that("residuals2 in mixed models", {

    ## mean
    expect_equal(colSums(e2.lvm$sCorrect$residuals),
                 colSums(residuals(e.lvm)))

    ## variance
    getVarCov2(e.lvm)/var(e2.lvm$sCorrect$residuals)
    ## not 1....
    
})

##----------------------------------------------------------------------
### test1-sCorrect-adjustedResiduals.R ends here

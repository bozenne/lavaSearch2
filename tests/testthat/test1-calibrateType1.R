### test1-calibrateType1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 10 2018 (09:34) 
## Version: 
## Last-Updated: apr 13 2018 (15:17) 
##           By: Brice Ozenne
##     Update #: 10
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
rm(list = ls())
if(FALSE){ ## already called in test-all.R    
    library(testthat)
    library(lavaSearch2)
}
lava.options(symbols = c("~","~~"))
context("calibrateType1")

## * Simulation
m.generative <- lvm(Y~X1+X2,G~1)
generative.coef <- c("Y" = 0.35,
                     "Y~X1" = 1.5,
                     "Y~X2" = 1.25,
                     "Y~~Y" = 1.14)
true.coef <- c("Y~G" = 0,
               generative.coef)
m.fit <- lvm(Y~X1+X2+G)

test_that("calibrateType1", {

    out <- calibrateType1(m.fit, true.coef = true.coef,
                          null = c("Y~G"), checkType1 = TRUE,
                          n = c(10,20), n.rep = 2,
                          generative.object = m.generative,
                          generative.coef = generative.coef,
                          dir.save = NULL, Ftest = TRUE,
                          bootstrap = FALSE, seed = 10, trace = 0)

    expect_equal(out$p.value[1:2,"p.Ztest"], rep(1.793176e-05,2))
    expect_equal(out$p.value[1:2,"p.robustZtest"], rep(0.0002798229,2))
    expect_equal(out$p.value[1:2,"p.Satt"], rep(0.001588207,2))
    expect_equal(out$p.value[1:2,"p.robustSatt"], rep(0.004587442,2))
    expect_equal(out$p.value[1:2,"p.SSC"], rep(0.0008923767,2))
    expect_equal(out$p.value[1:2,"p.robustSSC"], rep(0.004887758,2))
    expect_equal(out$p.value[1:2,"p.KR"], rep(0.01595612,2))
    expect_equal(out$p.value[1:2,"p.robustKR"], rep(0.03058102,2))
})


######################################################################
### test1-calibrateType1.R ends here

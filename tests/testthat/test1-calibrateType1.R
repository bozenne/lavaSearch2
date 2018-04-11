### test1-calibrateType1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 10 2018 (09:34) 
## Version: 
## Last-Updated: apr 11 2018 (13:07) 
##           By: Brice Ozenne
##     Update #: 9
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
                          dir.save = NULL,
                          bootstrap = FALSE, seed = 10, trace = 0)

    expect_equal(out$p.value[1,"p.Ztest"], 1.793176e-05)
    expect_equal(out$p.value[1,"p.robustZtest"], 0.0002798229)
    expect_equal(out$p.value[1,"p.Satt"], 0.001588207)
    expect_equal(out$p.value[1,"p.robustSatt"], 0.004587442)
    expect_equal(out$p.value[1,"p.SSC"], 0.0008923767)
    expect_equal(out$p.value[1,"p.robustSSC"], 0.004887758)
    expect_equal(out$p.value[1,"p.KR"], 0.01595612)
    expect_equal(out$p.value[1,"p.robustKR"], 0.03058102)
})


######################################################################
### test1-calibrateType1.R ends here

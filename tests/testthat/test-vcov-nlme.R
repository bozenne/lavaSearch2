### test-vcov-nlme.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec 13 2017 (15:52) 
## Version: 
## Last-Updated: dec 13 2017 (16:19) 
##           By: Brice Ozenne
##     Update #: 4
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

library(testthat)
library(nlme)

context("score2")

## * Corrected vcov
n <- 2e1

## ** linear model
m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e.gls <- gls(Y~X1+X2+X3, data = d, method = "ML")
e.lm <- lm(Y~X1+X2+X3, data = d)
Sigma.gls <- attr(residuals2(e.gls, cluster = 1:n, return.vcov.param = TRUE), "vcov.param")
Sigma.lm <- attr(residuals2(e.lm, return.vcov.param = TRUE), "vcov.param")


test_that("Corrected vcov - linear model",{
    expect_equal(unname(Sigma.gls), unname(Sigma.lm))
})


##----------------------------------------------------------------------
### test-vcov-nlme.R ends here

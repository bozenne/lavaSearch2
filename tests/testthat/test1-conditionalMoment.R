### test1-conditionalMoment.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 27 2018 (09:50) 
## Version: 
## Last-Updated: dec 17 2019 (13:47) 
##           By: Brice Ozenne
##     Update #: 25
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
context("conditionalMoment")

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

## * multiple linear regression
## ** no constrains
e.lvm <- estimate(lvm(Y1~X1,Y2~X2,Y3~X1+X3), data = d)

test_that("linear regression - no constrains",{
    test <- sCorrect(e.lvm, ssc = NULL, df = NULL)

    expect_equal(lava::score(e.lvm, indiv = TRUE),
                 score2(test, ssc = NULL, indiv = TRUE),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test, ssc = NULL)),
                 tol = 1e-8)
    expect_equal(coef(e.lvm),
                 coef2(test, ssc = NULL),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test, ssc = NULL, type = "response")),
                 tol = 1e-8)
})

## ** constrains and covariance link
m <- lvm(Y1[mu1:sigma]~X1,
         Y2[mu2:sigma]~X2,Y3~X1+X3,
         Y2~~Y3)
e.lvm <- estimate(m, data = d)

test_that("linear regression - constrains and covariance",{
    test <- sCorrect(e.lvm, ssc = NULL, df = NULL)
    
    expect_equal(unname(lava::score(e.lvm, indiv = TRUE)),
                 unname(score2(test, ssc = NULL, indiv = TRUE)),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test, ssc = NULL)),
                 tol = 1e-8)
    expect_equal(unname(coef(e.lvm)),
                 unname(coef2(test, ssc = NULL)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test, ssc = NULL, type = "response")),
                 tol = 1e-8)
})

## * factor model
m <- lvm(Y1~eta,
         Y2~eta+X2,
         Y3~eta,
         Z1~eta, Z1~~Y1,Z1~~Y2,
         eta~X1+X3)
e.lvm <- estimate(m, d)

test_that("factor model",{
    test <- sCorrect(e.lvm, ssc = NULL, df = NULL)
    
    expect_equal(lava::score(e.lvm, indiv = TRUE),
                 score2(test, ssc = NULL, indiv = TRUE),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test, ssc = NULL)),
                 tol = 1e-8)
    expect_equal(coef(e.lvm),
                 coef2(test, ssc = NULL),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test, ssc = NULL, type = "response")),
                 tol = 1e-8)
})

## * two factor model
## ** correlation
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,
           eta1~eta2))

e.lvm <- estimate(m, d)

test_that("two factor model - correlation",{
    test <- sCorrect(e.lvm, ssc = NULL, df = NULL)
    
    expect_equal(lava::score(e.lvm, indiv = TRUE),
                 score2(test, ssc = NULL, indiv = TRUE),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test, ssc = NULL)),
                 tol = 1e-8)
    expect_equal(coef(e.lvm),
                 coef2(test, ssc = NULL),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test, ssc = NULL, type = "response")),
                 tol = 1e-8)
})

## ** covariance
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,
           eta1~~eta2))

e.lvm <- estimate(m, d)

test_that("two factor model - covariance",{
    test <- sCorrect(e.lvm, ssc = NULL, df = NULL)
    
    expect_equal(lava::score(e.lvm, indiv = TRUE),
                 score2(test, ssc = NULL, indiv = TRUE),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test, ssc = NULL)),
                 tol = 1e-8)
    expect_equal(coef(e.lvm),
                 coef2(test, ssc = NULL),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test, ssc = NULL, type = "response")),
                 tol = 1e-8)
})

##
##----------------------------------------------------------------------
### test1-conditionalMoment.R ends here

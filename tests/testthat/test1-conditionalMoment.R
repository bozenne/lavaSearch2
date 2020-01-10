### test1-conditionalMoment.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 27 2018 (09:50) 
## Version: 
## Last-Updated: jan 10 2020 (14:12) 
##           By: Brice Ozenne
##     Update #: 42
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

test.secondOrder <- TRUE

## * simulation
n <- 5e1
mSim <- lvm(c(Y1~eta1,Y2~eta1+X2,Y3~eta1+X1,
              Z1~eta2,Z2~eta2,Z3~eta2+X3,
              X4~1,X5~1))
regression(mSim) <- eta1~X1+Gender
latent(mSim) <- ~eta1+eta2
categorical(mSim, labels = c("Male","Female")) <- ~Gender
transform(mSim, Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
d <- lava::sim(mSim, n = n, latent = FALSE)
dL <- reshape2::melt(d, id.vars = c("Id","X1","X2","X3","X4","X5","Gender"),
                     measure.vars = c("Y1","Y2","Y3","Z1","Z2","Z3"))

## * multiple linear regression
    
## ** no constrains
e.lvm <- estimate(lvm(Y1~X1,Y2~X2,Y3~X1+X3), data = d)
e.lm <- lm(Y1~X1, data = d)
lvm2lm <- c("(Intercept)" = "Y1", "X1" = "Y1~X1", "sigma2" = "Y1~~Y1")

test_that("linear regression - no constrains",{
    test.lvm <- sCorrect(e.lvm, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")
    test.lm <- sCorrect(e.lm, param = setNames(coef(e.lvm)[lvm2lm], names(lvm2lm)),
                        ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")
    
    expect_equal(lava::score(e.lvm, indiv = TRUE),
                 score2(test.lvm, ssc = NA, indiv = TRUE),
                 tol = 1e-8)
    expect_equal(unname(score2(test.lm, ssc = NA, indiv = TRUE)),
                 unname(score2(test.lvm, ssc = NA, indiv = TRUE)[,lvm2lm]),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(information(test.lm, ssc = NA)),
                 unname(information2(test.lvm, ssc = NA)[lvm2lm,lvm2lm]),
                 tol = 1e-8)
    expect_equal(coef(e.lvm),
                 coef2(test.lvm, ssc = NA),
                 tol = 1e-8)
    expect_equal(unname(coef2(test.lm, ssc = NA)[1:3]),
                 unname(coef2(test.lvm, ssc = NA)[lvm2lm]),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)
    expect_equal(unname(residuals2(test.lm, ssc = NA, type = "response")[,1]),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")[,"Y1"]),
                 tol = 1e-8)

    if(test.secondOrder){
        testN.lvm <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "numeric")
        expect_equal(test.lvm$sCorrect$hessian,testN.lvm$sCorrect$hessian, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dVcov.param,
                     testN.lvm$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dRvcov.param,testN.lvm$sCorrect$dRvcov.param, tol = 1e-7)

        expect_equal(unname(test.lvm$sCorrect$hessian[lvm2lm,lvm2lm,]),unname(test.lm$sCorrect$hessian), tol = 1e-7)
        expect_equal(unname(test.lvm$sCorrect$dVcov.param[lvm2lm,lvm2lm,lvm2lm,drop=FALSE]),
                     unname(test.lm$sCorrect$dVcov.param), tol = 1e-7)
        expect_equal(unname(test.lvm$sCorrect$dRvcov.param[lvm2lm,lvm2lm,lvm2lm,drop=FALSE]),
                     unname(test.lm$sCorrect$dRvcov.param), tol = 1e-7)
    }
})

## ** constrains and covariance link
m <- lvm(Y1[mu1:sigma]~X1,
         Y2[mu2:sigma]~X2,
         Y3~X1+X3,
         Y2~~Y3)
e.lvm <- estimate(m, data = d)

test_that("linear regression - constrains and covariance",{
    test <- sCorrect(e.lvm, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")
    
    expect_equal(unname(lava::score(e.lvm, indiv = TRUE)),
                 unname(score2(test, ssc = NA, indiv = TRUE)),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(coef(e.lvm)),
                 unname(coef2(test, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test, ssc = NA, type = "response")),
                 tol = 1e-8)

    if(test.secondOrder){
        testN <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "numeric")
        keep.param <- dimnames(test$sCorrect$dVcov.param)[[3]]
        zero.param <- setdiff(dimnames(testN$sCorrect$dVcov.param)[[3]], keep.param)
        expect_equal(test$sCorrect$hessian,testN$sCorrect$hessian, tol = 1e-7)
        expect_equal(test$sCorrect$dVcov.param,testN$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test$sCorrect$dRvcov.param,testN$sCorrect$dRvcov.param, tol = 1e-7)
    }
})

## ** gls with heterogenous variance
e.gls <- gls(value ~ 0+variable*X1 + X2, weight = varIdent(form =~ 1 | variable), data = dL[dL$variable %in% c("Y1","Y2"),],
             method = "ML")
e.lvm <- estimate(lvm(Y1[mu1:sigma1] ~ X1 + b*X2, Y2[mu2:sigma2] ~ X1 + b*X2), data = d)
## logLik(e.gls)
## logLik(e.lvm)

test_that("gls",{
    test.gls <- sCorrect(e.gls, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")
    test.lvm <- sCorrect(e.lvm, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")

    expect_equal(unname(lava::score(e.lvm, indiv = TRUE)),
                 unname(score2(test.lvm, ssc = NA, indiv = TRUE)),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(coef(e.lvm)),
                 unname(coef2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)

    if(test.secondOrder){
        testN.lvm <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "numeric")

        expect_equal(test.lvm$sCorrect$hessian,testN.lvm$sCorrect$hessian, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dVcov.param,testN.lvm$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dRvcov.param,testN.lvm$sCorrect$dRvcov.param, tol = 1e-7)
    }
})


## * factor model
m <- lvm(Y1~eta,
         Y2~eta+X2,
         Y3~eta,
         Z1~eta, Z1~~Y1,Z1~~Y2,
         eta~X1+X3)
e.lvm <- estimate(m, d)

test_that("factor model",{
    test <- sCorrect(e.lvm, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")
    
    expect_equal(lava::score(e.lvm, indiv = TRUE),
                 score2(test, ssc = NA, indiv = TRUE),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test, ssc = NA)),
                 tol = 1e-8)
    expect_equal(coef(e.lvm),
                 coef2(test, ssc = NA),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test, ssc = NA, type = "response")),
                 tol = 1e-8)

    if(test.secondOrder){
        testN <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "numeric")

        expect_equal(test$sCorrect$hessian,testN$sCorrect$hessian, tol = 1e-7)
        expect_equal(test$sCorrect$dVcov.param,testN$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test$sCorrect$dRvcov.param,testN$sCorrect$dRvcov.param, tol = 1e-7)
    }
})

## * two factor model
## ** correlation
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,
           eta1~eta2))

e.lvm <- estimate(m, d)

test_that("two factor model - correlation",{
    test <- sCorrect(e.lvm, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")
    
    expect_equal(lava::score(e.lvm, indiv = TRUE),
                 score2(test, ssc = NA, indiv = TRUE),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test, ssc = NA)),
                 tol = 1e-8)
    expect_equal(coef(e.lvm),
                 coef2(test, ssc = NA),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test, ssc = NA, type = "response")),
                 tol = 1e-8)

    if(test.secondOrder){
        testN <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "numeric")

        expect_equal(test$sCorrect$hessian,testN$sCorrect$hessian, tol = 1e-7)
        expect_equal(test$sCorrect$dVcov.param,testN$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test$sCorrect$dRvcov.param,testN$sCorrect$dRvcov.param, tol = 1e-7)
    }
})

## ** covariance
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,eta1~X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,eta2~X4,
           eta1~~eta2))

e.lvm <- estimate(m, d)

test_that("two factor model - covariance",{
    test <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "analytic")
    
    expect_equal(lava::score(e.lvm, indiv = TRUE),
                 score2(test, ssc = NA, indiv = TRUE),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm),
                 unname(information2(test, ssc = NA)),
                 tol = 1e-8)
    expect_equal(coef(e.lvm),
                 coef2(test, ssc = NA),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm)),
                 unname(residuals2(test, ssc = NA, type = "response")),
                 tol = 1e-8)

    if(test.secondOrder){
        testN <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "numeric")
    
        expect_equal(test$sCorrect$hessian,testN$sCorrect$hessian, tol = 1e-7)
        expect_equal(test$sCorrect$dVcov.param,testN$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test$sCorrect$dRvcov.param,testN$sCorrect$dRvcov.param, tol = 1e-7)
    }

})

##
##----------------------------------------------------------------------
### test1-conditionalMoment.R ends here

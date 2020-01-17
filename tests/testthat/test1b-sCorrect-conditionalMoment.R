### test1-conditionalMoment.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar 27 2018 (09:50) 
## Version: 
## Last-Updated: jan 16 2020 (15:57) 
##           By: Brice Ozenne
##     Update #: 71
##----------------------------------------------------------------------
## 
### Commentary: 
## Compare the computation of the score/information matrix/residuals between sCorrect and lava
## Compare the computation of the hessian/derivative of the information matrix using analytical formulae vs. numerical derivatives
##
## NOTE: iid in lava uses numerical derivative to compute the information matrix
## this is why there is not a perfect matching between iid2.lvm and iid.lvm
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
a
## * simulation
cat("- simulation \n")
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
cat("- multiple linear regression \n")
    
## ** no constrains
e.lvm <- estimate(lvm(Y1~X1,Y2~X2,Y3~X1+X3), data = d)
e.lm <- lm(Y1~X1, data = d)
lvm2lm <- c("(Intercept)" = "Y1", "X1" = "Y1~X1", "sigma2" = "Y1~~Y1")

test_that("linear regression (ML) - no constrains",{
    expect_equal(as.double(logLik(e.lvm)), -277.9339, tol = 1e-6)
    
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
        
        X <- lapply(list(Y1~X1,Y2~X2,Y3~X1+X3), model.matrix, d)
        sigma2 <- list(coef(e.lvm)["Y1~~Y1"],
                       coef(e.lvm)["Y2~~Y2"],
                       coef(e.lvm)["Y3~~Y3"])

        dI <- mapply(X,sigma2, FUN = function(iX,iSigma){
            bdiag(crossprod(iX)/iSigma^2,n/(iSigma^3))
        })
        vcov <- mapply(X,sigma2, FUN = function(iX,iSigma){
            solve(bdiag(crossprod(iX)/iSigma,n/(2*iSigma^2)))
        })
        GS <- mapply(vcov, dI, FUN = function(x,y){
            as.matrix(x %*% y %*% x)
        })
        nameCoef.Y <- lapply(list("Y1","Y2","Y3"), function(iY){grep(iY, names(coef(e.lvm)), value = TRUE)})

        testN.lvm <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "numeric")
        expect_equal(test.lvm$sCorrect$hessian,testN.lvm$sCorrect$hessian, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dVcov.param, testN.lvm$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(unname(test.lvm$sCorrect$dVcov.param[nameCoef.Y[[1]],nameCoef.Y[[1]],"Y1~~Y1"]), unname(GS[[1]]), tol = 1e-7)
        expect_equal(unname(test.lvm$sCorrect$dVcov.param[nameCoef.Y[[2]],nameCoef.Y[[2]],"Y2~~Y2"]), unname(GS[[2]]), tol = 1e-7)
        expect_equal(unname(test.lvm$sCorrect$dVcov.param[nameCoef.Y[[3]],nameCoef.Y[[3]],"Y3~~Y3"]), unname(GS[[3]]), tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dRvcov.param,testN.lvm$sCorrect$dRvcov.param, tol = 1e-7)

        expect_equal(unname(test.lvm$sCorrect$hessian[lvm2lm,lvm2lm,]),unname(test.lm$sCorrect$hessian), tol = 1e-7)
        expect_equal(unname(test.lvm$sCorrect$dVcov.param[lvm2lm,lvm2lm,lvm2lm,drop=FALSE]),
                     unname(test.lm$sCorrect$dVcov.param), tol = 1e-7)
        expect_equal(unname(test.lvm$sCorrect$dRvcov.param[lvm2lm,lvm2lm,lvm2lm,drop=FALSE]),
                     unname(test.lm$sCorrect$dRvcov.param), tol = 1e-7)
    }
})

test_that("linear regression (ML+1) - no constrains",{
    newcoef <- coef(e.lvm)+1
    test.lvm <- sCorrect(e.lvm, param = newcoef, ssc = NA, df = NA)
    
    expect_equal(lava::score(e.lvm, p = newcoef, indiv = TRUE),
                 score2(test.lvm, ssc = NA, indiv = TRUE),
                 tol = 1e-8)
    expect_equal(lava::information(e.lvm, p = newcoef),
                 unname(information2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm, p = newcoef)),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)
})

## ** constrains and covariance link
m <- lvm(Y1[mu1:sigma]~X1,
         Y2[mu2:sigma]~X2,
         Y3~X1+X3,
         Y2~~Y3)
e.lvm <- estimate(m, data = d)

test_that("linear regression - constrains and covariance",{
    expect_equal(as.double(logLik(e.lvm)), -272.5088, tol = 1e-6)

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

        ## compare to previous versions
        GS <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.00118648, 0.01999964, 0, -0.00554516, 3.94e-06, -3.01e-05, 0, 0, 0, 0, 0.01999964, -4.014e-05, 0, -1.68e-06, 0.00086362, -3.14e-05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.00554516, -1.68e-06, 0, -0.02591598, 1.841e-05, -0.00014069, 0, 0, 0, 0, 3.94e-06, 0.00086362, 0, 1.841e-05, -0.01884834, -0.00143231, 0, 0, 0, 0, -3.01e-05, -3.14e-05, 0, -0.00014069, -0.00143231, -0.01661326, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05139595, 0.05281454, 0, 0, 0, 0, 0, 0, 0, 0.05139595, -0.02433603, 0.0498916, 0, 0, 0, 0, 0, 0, 0, 0.05281454, 0.0498916, 0), 
                     nrow = 10, 
                     ncol = 10, 
                     dimnames = list(c("mu1", "mu2", "Y3", "Y1~X1", "Y2~X2", "Y3~X1", "Y3~X3", "sigma", "Y3~~Y3", "Y2~~Y3"),c("mu1", "mu2", "Y3", "Y1~X1", "Y2~X2", "Y3~X1", "Y3~X3", "sigma", "Y3~~Y3", "Y2~~Y3")) 
                     ) 
        expect_equal(test$sCorrect$dVcov.param[,,"Y2~~Y3"],GS,tol = 1e-6)
    }
})

test_that("linear regression (ML+1) - constrains and covariance",{
    newcoef <- coef(e.lvm)+1
    test.lvm <- sCorrect(e.lvm, param = .coef2(e.lvm, labels = 1, ssc = NA)+1, ssc = NA, df = NA)
    
    expect_equal(unname(lava::score(e.lvm, p = newcoef, indiv = TRUE)),
                 unname(score2(test.lvm, ssc = NA, indiv = TRUE)),
                 tol = 1e-8)
    expect_equal(unname(lava::information(e.lvm, p = newcoef)),
                 unname(information2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm, p = newcoef)),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)
})

## ** gls with heterogenous variance
e.gls <- gls(value ~ 0+variable*X1 + X2, weight = varIdent(form =~ 1 | variable), data = dL[dL$variable %in% c("Y1","Y2"),],
             method = "ML")
e.lvm <- estimate(lvm(Y1[mu1:sigma1] ~ X1 + b*X2, Y2[mu2:sigma2] ~ X1 + b*X2), data = d)
## logLik(e.gls)
## logLik(e.lvm)

test_that("gls",{
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.gls)), tol = 1e-3)
    expect_equal(as.double(logLik(e.lvm)), -179.2554, tol = 1e-6)

    test.gls <- sCorrect(e.gls, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")
    test.lvm <- sCorrect(e.lvm, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")

    expect_equal(unname(getVarCov2(test.gls, ssc = NA)),
                 unname(getVarCov2(test.lvm, ssc = NA)),
                 tol = 1e-3)
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
    test <- residuals2(test.gls, ssc = NA, type = "response")
    GS <- residuals2(test.lvm, ssc = NA, type = "response")
    expect_equal(as.vector(na.omit(as.double(test))),
                 as.double(GS),
                 tol = 1e-4)
    
    if(test.secondOrder){
        testN.gls <- sCorrect(e.gls, ssc = NA, df = "Satterthwaite", derivative = "numeric")
        testN.lvm <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "numeric")

        expect_equal(test.lvm$sCorrect$hessian,testN.lvm$sCorrect$hessian, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dVcov.param,testN.lvm$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dRvcov.param,testN.lvm$sCorrect$dRvcov.param, tol = 1e-7)

        expect_equal(test.gls$sCorrect$hessian,testN.gls$sCorrect$hessian, tol = 1e-7)
        expect_equal(test.gls$sCorrect$dVcov.param,testN.gls$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test.gls$sCorrect$dRvcov.param,testN.gls$sCorrect$dRvcov.param, tol = 1e-7)

        ## compare to previous versions
        GS <- matrix(c(0.00027028, 0.00027028, -1.22e-05, 0.00126055, -1.22e-05, 0, 0, 0.00027028, 0.02031147, -1.22e-05, 0.00126055, -0.00090249, 0, 0, -1.22e-05, -1.22e-05, 5.5e-07, -5.691e-05, 5.5e-07, 0, 0, 0.00126055, 0.00126055, -5.691e-05, 0.00587905, -5.691e-05, 0, 0, -1.22e-05, -0.00090249, 5.5e-07, -5.691e-05, 0.01924312, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.16350761), 
                     nrow = 7, 
                     ncol = 7, 
                     dimnames = list(c("mu1", "mu2", "Y1~X1", "b", "Y2~X1", "sigma1", "sigma2"),c("mu1", "mu2", "Y1~X1", "b", "Y2~X1", "sigma1", "sigma2")) 
                     ) 
        expect_equal(test.lvm$sCorrect$dVcov.param[,,"sigma2"],GS,tol = 1e-6)
    }
})

test_that("linear regression (ML+1) - gls",{
    newcoef <- coef(e.lvm)+1
    test.lvm <- sCorrect(e.lvm, param = .coef2(e.lvm, labels = 1, ssc = NA)+1, ssc = NA, df = NA)
    
    expect_equal(unname(lava::score(e.lvm, p = newcoef, indiv = TRUE)),
                 unname(score2(test.lvm, ssc = NA, indiv = TRUE)),
                 tol = 1e-8)
    expect_equal(unname(lava::information(e.lvm, p = newcoef)),
                 unname(information2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm, p = newcoef)),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)
})

## * mixed model
cat("- mixed model \n")

## ** Compound symmetry
m <- lvm(Y1[mu1:sigma]~1*eta,
         Y2[mu2:sigma]~1*eta,
         Y3[mu3:sigma]~1*eta,
         eta~X1+Gender)
e.lvm <- estimate(m, d)

e.lme <- lme(value ~ variable + X1 + Gender,
             random =~ 1|Id,
             data = dL[dL$variable %in% c("Y1","Y2","Y3"),],
             method = "ML")

e.gls <- gls(value ~ variable + X1 + Gender,
             correlation = corCompSymm(form=~ 1|Id),
             data = dL[dL$variable %in% c("Y1","Y2","Y3"),],
             method = "ML")

test_that("Compound symmetry", {
    expect_equal(as.double(logLik(e.lme)),as.double(logLik(e.lvm)), tol = 1e-3)
    expect_equal(as.double(logLik(e.gls)),as.double(logLik(e.lvm)), tol = 1e-3)
    expect_equal(as.double(logLik(e.lvm)), -259.8317, tol = 1e-6)

    test.gls <- sCorrect(e.gls, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")
    test.lme <- sCorrect(e.lme, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")
    test.lvm <- sCorrect(e.lvm, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")

    expect_equal(unname(getVarCov2(test.gls, ssc = NA)),
                 unname(getVarCov2(test.lvm, ssc = NA)),
                 tol = 1e-3)
    expect_equal(unname(getVarCov2(test.lme, ssc = NA)),
                 unname(getVarCov2(test.lvm, ssc = NA)),
                 tol = 1e-3)
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
    expect_equal(unname(residuals2(test.gls, ssc = NA, type = "response")),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)

    if(test.secondOrder){
        testN.gls <- sCorrect(e.gls, ssc = NA, df = "Satterthwaite", derivative = "numeric")
        testN.lvm <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "numeric")

        expect_equal(test.gls$sCorrect$hessian,testN.gls$sCorrect$hessian, tol = 1e-7)
        expect_equal(test.gls$sCorrect$dVcov.param,testN.gls$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test.gls$sCorrect$dRvcov.param,testN.gls$sCorrect$dRvcov.param, tol = 1e-7)

        expect_equal(test.lvm$sCorrect$hessian,testN.lvm$sCorrect$hessian, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dVcov.param,testN.lvm$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dRvcov.param,testN.lvm$sCorrect$dRvcov.param, tol = 1e-7)

        ## compare to previous versions
        GS <- matrix(c(0.02860806, -0.02, -0.02, -0.00089185, -0.01529786, 0, 0, -0.02, 0.04, 0.02, 0, 0, 0, 0, -0.02, 0.02, 0.04, 0, 0, 0, 0, -0.00089185, 0, 0, 0.00645539, 0.00105925, 0, 0, -0.01529786, 0, 0, 0.00105925, 0.02723009, 0, 0, 0, 0, 0, 0, 0, 0.05700596, -0.01900199, 0, 0, 0, 0, 0, -0.01900199, 0.03500504), 
                     nrow = 7, 
                     ncol = 7, 
                     dimnames = list(c("eta", "mu2", "mu3", "eta~X1", "eta~GenderFemale", "sigma", "eta~~eta"),c("eta", "mu2", "mu3", "eta~X1", "eta~GenderFemale", "sigma", "eta~~eta")) 
                     )
        expect_equal(test.lvm$sCorrect$dVcov.param[,,"sigma"],GS,tol = 1e-6)
    }
})

test_that("mixed model (ML+1) - CS",{
    newcoef <- coef(e.lvm)+1
    test.lvm <- sCorrect(e.lvm, param = .coef2(e.lvm, labels = 1, ssc = NA)+1, ssc = NA, df = NA)
    
    expect_equal(unname(lava::score(e.lvm, p = newcoef, indiv = TRUE)),
                 unname(score2(test.lvm, ssc = NA, indiv = TRUE)),
                 tol = 1e-8)
    expect_equal(unname(lava::information(e.lvm, p = newcoef)),
                 unname(information2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm, p = newcoef)),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)
})

## ** Unstructured 
m <- lvm(Y1~1*eta,
         Y2~1*eta,
         Y3~1*eta,
         eta~X1+Gender)
covariance(m) <- Y1~Y2
covariance(m) <- Y1~Y3
e.lvm <- estimate(m, d)

e.lme <- lme(value ~ variable + X1 + Gender,
             random =~ 1|Id,
             correlation = corSymm(),
             weights = varIdent(form =~ 1|variable),
             data = dL[dL$variable %in% c("Y1","Y2","Y3"),],
             method = "ML")

e.gls <- gls(value ~ variable + X1 + Gender,
             correlation = corSymm(form=~ 1|Id),
             weights = varIdent(form =~ 1|variable),
             data = dL[dL$variable %in% c("Y1","Y2","Y3"),],
             method = "ML")

test_that("Unstructured", {
    expect_equal(as.double(logLik(e.lme)),as.double(logLik(e.lvm)), tol = 1e-3)
    expect_equal(as.double(logLik(e.gls)),as.double(logLik(e.lvm)), tol = 1e-3)
    expect_equal(as.double(logLik(e.lvm)), -258.8121, tol = 1e-6)

    test.gls <- sCorrect(e.gls, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")
    test.lme <- sCorrect(e.lme, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic") ## error: overparametrized model
    test.lvm <- sCorrect(e.lvm, ssc = NA, df = if(test.secondOrder){"Satterthwaite"}else{NA}, derivative = "analytic")

    expect_equal(unname(getVarCov2(test.gls, ssc = NA)),
                 unname(getVarCov2(test.lvm, ssc = NA)),
                 tol = 1e-3)
    expect_equal(unname(getVarCov2(test.lme, ssc = NA)),
                 unname(getVarCov2(test.lvm, ssc = NA)),
                 tol = 1e-3)
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
    expect_equal(unname(residuals2(test.gls, ssc = NA, type = "response")),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-3)

    if(test.secondOrder){
        testN.gls <- sCorrect(e.gls, ssc = NA, df = "Satterthwaite", derivative = "numeric")
        testN.lvm <- sCorrect(e.lvm, ssc = NA, df = "Satterthwaite", derivative = "numeric")

        expect_equal(test.gls$sCorrect$hessian,testN.gls$sCorrect$hessian, tol = 1e-7)
        expect_equal(test.gls$sCorrect$dVcov.param,testN.gls$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test.gls$sCorrect$dRvcov.param,testN.gls$sCorrect$dRvcov.param, tol = 1e-7)

        expect_equal(test.lvm$sCorrect$hessian,testN.lvm$sCorrect$hessian, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dVcov.param,testN.lvm$sCorrect$dVcov.param, tol = 1e-7)
        expect_equal(test.lvm$sCorrect$dRvcov.param,testN.lvm$sCorrect$dRvcov.param, tol = 1e-7)

        ## compare to previous versions
        GS <- matrix(c(0.02337325, -0.02, -0.02, -0.00034949, -0.00599479, 0, 0, 0, 0, 0, 0, -0.02, 0.02, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, -0.02, 0.02, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, -0.00034949, 0, 0, 0.00252968, 0.00041509, 0, 0, 0, 0, 0, 0, -0.00599479, 0, 0, 0.00041509, 0.0106707, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.14197442, 0, 0, 0, 0.02650713, 0.02772787, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.02650713, 0, 0, 0, 0.04155753, 0.00895138, 0, 0, 0, 0, 0, 0.02772787, 0, 0, 0, 0.00895138, 0.04452667), 
                     nrow = 11, 
                     ncol = 11, 
                     dimnames = list(c("eta", "Y2", "Y3", "eta~X1", "eta~GenderFemale", "Y1~~Y1", "eta~~eta", "Y2~~Y2", "Y3~~Y3", "Y1~~Y2", "Y1~~Y3"),c("eta", "Y2", "Y3", "eta~X1", "eta~GenderFemale", "Y1~~Y1", "eta~~eta", "Y2~~Y2", "Y3~~Y3", "Y1~~Y2", "Y1~~Y3")) 
                     ) 
        expect_equal(test.lvm$sCorrect$dVcov.param[,,"Y1~~Y1"],GS,tol = 1e-6)
    }

})

test_that("mixed model (ML+1) - UN",{
    newcoef <- coef(e.lvm)+1
    test.lvm <- sCorrect(e.lvm, param = .coef2(e.lvm, labels = 1, ssc = NA)+1, ssc = NA, df = NA)
    
    expect_equal(unname(lava::score(e.lvm, p = newcoef, indiv = TRUE)),
                 unname(score2(test.lvm, ssc = NA, indiv = TRUE)),
                 tol = 1e-8)
    expect_equal(unname(lava::information(e.lvm, p = newcoef)),
                 unname(information2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm, p = newcoef)),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)
})

## * factor model
cat("- factor model \n")

m <- lvm(Y1~eta,
         Y2~eta+X2,
         Y3~eta,
         Z1~eta, Z1~~Y1,Z1~~Y2,
         eta~X1+X3)
e.lvm <- estimate(m, d)

test_that("factor model",{
    expect_equal(as.double(logLik(e.lvm)), -334.2583, tol = 1e-6)

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

        ## compare to previous versions
        GS <- matrix(c(0.04792556, 0.04673097, 0.04673097, -0.00085306), 
                     nrow = 2, 
                     ncol = 2, 
                     dimnames = list(c("Y2", "Y3"),c("Y2", "Y3")) 
                     ) 
        expect_equal(test$sCorrect$dVcov.param[paste0("Y",2:3),paste0("Y",2:3),"Y2~eta"],GS,tol = 1e-6)
    }
})

test_that("factor model (ML+1)",{
    newcoef <- coef(e.lvm)+1
    test.lvm <- sCorrect(e.lvm, param = .coef2(e.lvm, labels = 1, ssc = NA)+1, ssc = NA, df = NA)
    
    expect_equal(unname(lava::score(e.lvm, p = newcoef, indiv = TRUE)),
                 unname(score2(test.lvm, ssc = NA, indiv = TRUE)),
                 tol = 1e-8)
    expect_equal(unname(lava::information(e.lvm, p = newcoef)),
                 unname(information2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm, p = newcoef)),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)
})

## * two factor model
cat("- two factor model \n")

## ** correlation
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,
           eta1~eta2))

e.lvm <- estimate(m, d)

test_that("two factor model - correlation",{
    expect_equal(as.double(logLik(e.lvm)), -518.3538, tol = 1e-6)

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
    
        ## compare to previous versions
        GS <- matrix(c(-0.01772622, 0.01772622, -0.01540713, 0.07276457, 0.01772622, -0.01772622, 0.01540713, -0.07276457, -0.01540713, 0.01540713, -0.0133915, 0.06324503, 0.07276457, -0.07276457, 0.06324503, -0.29869238), 
                     nrow = 4, 
                     ncol = 4, 
                     dimnames = list(c("Z1~~Z1", "eta2~~eta2", "Z2~~Z2", "Z3~~Z3"),c("Z1~~Z1", "eta2~~eta2", "Z2~~Z2", "Z3~~Z3")) 
                     )
        keep.param <-  c("Z1~~Z1", "eta2~~eta2", "Z2~~Z2", "Z3~~Z3")
        expect_equal(test$sCorrect$dVcov.param[keep.param,keep.param,"eta1~eta2"],GS,tol = 1e-6)
    }
})

test_that("two factor model (ML+1) - correlation",{
    newcoef <- coef(e.lvm)+1
    test.lvm <- sCorrect(e.lvm, param = .coef2(e.lvm, labels = 1, ssc = NA)+1, ssc = NA, df = NA)
    
    expect_equal(unname(lava::score(e.lvm, p = newcoef, indiv = TRUE)),
                 unname(score2(test.lvm, ssc = NA, indiv = TRUE)),
                 tol = 1e-8)
    expect_equal(unname(lava::information(e.lvm, p = newcoef)),
                 unname(information2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm, p = newcoef)),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)
})

## ** covariance
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,eta1~X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,eta2~X4,
           eta1~~eta2))

e.lvm <- estimate(m, d)

test_that("two factor model - covariance",{
    expect_equal(as.double(logLik(e.lvm)), -502.7317, tol = 1e-6)

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

        GS <- matrix(c(-0.01561465, 0.01447463, -0.01301386, 0.05822372, 0.01447463, -0.01341929, 0.01206254, -0.05397059, -0.01301386, 0.01206254, -0.01091509, 0.04865601, 0.05822372, -0.05397059, 0.04865601, -0.21734976), 
                     nrow = 4, 
                     ncol = 4, 
                     dimnames = list(c("Z1~~Z1", "eta2~~eta2", "Z2~~Z2", "Z3~~Z3"),c("Z1~~Z1", "eta2~~eta2", "Z2~~Z2", "Z3~~Z3")) 
                     ) 
        keep.param <-  c("Z1~~Z1", "eta2~~eta2", "Z2~~Z2", "Z3~~Z3")
        expect_equal(test$sCorrect$dVcov.param[keep.param,keep.param,"eta1~~eta2"],GS,tol = 1e-6)
    }
})

test_that("two factor model (ML+1) - covariance",{
    newcoef <- coef(e.lvm)+1
    test.lvm <- sCorrect(e.lvm, param = .coef2(e.lvm, labels = 1, ssc = NA)+1, ssc = NA, df = NA)
    
    expect_equal(unname(lava::score(e.lvm, p = newcoef, indiv = TRUE)),
                 unname(score2(test.lvm, ssc = NA, indiv = TRUE)),
                 tol = 1e-8)
    expect_equal(unname(lava::information(e.lvm, p = newcoef)),
                 unname(information2(test.lvm, ssc = NA)),
                 tol = 1e-8)
    expect_equal(unname(residuals(e.lvm, p = newcoef)),
                 unname(residuals2(test.lvm, ssc = NA, type = "response")),
                 tol = 1e-8)
})

##
##----------------------------------------------------------------------
### test1-conditionalMoment.R ends here

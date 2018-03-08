### test1-sCorrect-smallSampleCorrection.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  7 2018 (12:21) 
## Version: 
## Last-Updated: mar  8 2018 (16:37) 
##           By: Brice Ozenne
##     Update #: 7
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
if(FALSE){ ## already called in test-all.R
    rm(list = ls())
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))
context("sCorrect: small sample correction")

## * simulation
n <- 5e1
mSim <- lvm(c(Y1~eta1,Y2~eta1+X2,Y3~eta1+X1,
              Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(mSim) <- eta1~X1+Gender
latent(mSim) <- ~eta1+eta2
categorical(mSim, labels = c("Male","Female")) <- ~Gender
transform(mSim, Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
d <- sim(mSim, n = n, latent = FALSE)
dL <- melt(d, id.vars = c("Id","X1","X2","X3","Gender"),
           measure.vars = c("Y1","Y2","Y3","Z1","Z2","Z3"))
dLred <- dL[dL$variable %in% c("Y1","Y2","Y3"),]

## * linear regression [lm,gls,lvm]
## ** model fit and sCorrect
e.lvm <- estimate(lvm(Y1~X1+X2+Gender), data = d)
e.lm <- lm(Y1~X1+X2+Gender, data = d)
e.gls <- gls(Y1~X1+X2+Gender, data = d, method = "ML")

e2.lvm <- e.lvm
e2.gls <- e.gls
e2.lm <- e.lm



sCorrect(e2.lvm) <- FALSE
sCorrect(e2.gls, cluster = 1:n) <- FALSE
sCorrect(e2.lm) <- FALSE


## * adjusted residuals

## ** univariate linear model
m <- lvm(Y~X)
n <- 1e2
d <- sim(m,n)

test_that("residuals2 match residuals.lm (lm adjusted)", {
    e.lm <- lm(Y~X, data = d)
    epsilon.lm <- residuals(e.lm)
    X <- model.matrix(e.lm, d)
    iH <- diag(1,n,n) - X %*% solve(t(X) %*% X) %*% t(X)
    GS1 <- epsilon.lm/diag(iH)^(1/2)
    
    e.lvm <- estimate(lvm(Y~X), d)
    res2 <- residuals2(e.lvm, bias.correct = TRUE)

    expect_equal(as.double(res2),as.double(GS1))
})



## ** multivariate linear models
m <- lvm(Y~G+X,G~X)
n <- 1e2
d <- sim(m,n)

test_that("residuals2 match residuals.lm", {
    ## first model
    e.lm1 <- lm(Y~G+X, data = d)
    res1 <- residuals2(e.lm1, bias.correct = TRUE)

    epsilon.lm1 <- residuals(e.lm1)
    X1 <- model.matrix(e.lm1, d)
    iH1 <- diag(1,n,n) - X1 %*% solve(t(X1) %*% X1) %*% t(X1)

    expect_equal(as.double(epsilon.lm1/diag(iH1)^(1/2)),
                 as.double(res1))

    ## second model
    e.lm2 <- lm(G~X, data = d)
    res2 <- residuals2(e.lm2, bias.correct = TRUE)
    epsilon.lm2 <- residuals(e.lm2)
    X2 <- model.matrix(e.lm2, d)
    iH2 <- diag(1,n,n) - X2 %*% solve(t(X2) %*% X2) %*% t(X2)
    
    expect_equal(as.double(epsilon.lm2/diag(iH2)^(1/2)),
                 as.double(res2))

    ## global
    e.lvm <- estimate(m, d)    
    resTest <- residuals2(e.lvm, bias.correct = TRUE)

    expect_equal(as.double(resTest[,1]),as.double(res1))
    expect_equal(as.double(resTest[,2]),as.double(res2))   
})


## * Corrected vcov
n <- 2e1

## ** linear model
m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e.lvm <- estimate(lvm(Y~X1+X2+X3),d)
e.lm <- lm(Y~X1+X2+X3, data = d)
Sigma.lvm <- attr(residuals2(e.lvm, bias.correct = TRUE, return.vcov.param = TRUE), "vcov.param")
Sigma.lm <- attr(residuals2(e.lm, bias.correct = TRUE, return.vcov.param = TRUE), "vcov.param")

test_that("Corrected vcov - linear model",{
    expect_equal(unname(Sigma.lvm), unname(Sigma.lm))

    Sigma.uncorrected <- vcov(e.lvm)
    attr(Sigma.uncorrected, "det") <- NULL
    attr(Sigma.uncorrected, "pseudo") <- NULL
    attr(Sigma.uncorrected, "minSV") <- NULL

    p <- 4
    factor_beta <- (n+p)/n 
    expect_equal(Sigma.lvm[1:p,1:p], factor_beta * Sigma.uncorrected[1:p,1:p])

    factor_sigma <- n/(n-p)*(n+p)^2/n^2
    expect_equal(Sigma.lvm[p+1,p+1], factor_sigma * Sigma.uncorrected[p+1,p+1])
})

## ** multiple linear model

## *** simulation
m <- lvm(c(Y1,Y2)~X1+X2)
set.seed(10)
df <- sim(m,2e1)

## *** model fit
e1.lm <- lm(Y1~X1+X2, data = df)
e2.lm <- lm(Y2~X1+X2, data = df)
e.lvm <- estimate(m, data = df)

## *** tests
Sigma.lvm <- attr(residuals2(e.lvm, return.vcov.param = TRUE, bias.correct = TRUE),
                 "vcov.param")

Sigma.lvm[c("Y1","Y1~X1","Y1~X2"),c("Y1","Y1~X1","Y1~X2")]/vcov(e1.lm)
Sigma.lvm[c("Y2","Y2~X1","Y2~X2"),c("Y2","Y2~X1","Y2~X2")]/vcov(e2.lm)

Sigma.lvm[c("Y1~~Y1"),c("Y1~~Y1")]/(2*sigma(e1.lm)^4/(NROW(df)-length(coef(e1.lm))))
Sigma.lvm[c("Y2~~Y2"),c("Y2~~Y2")]/(2*sigma(e2.lm)^4/(NROW(df)-length(coef(e2.lm))))



## ** mixed model

## *** simulate
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G+Gender,X1~1,X2~1))
latent(mSim) <- ~eta
categorical(mSim, labels = c("M","F")) <- ~Gender
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- sim(mSim,n,latent = FALSE)
dW <- dW[order(dW$Id),,drop=FALSE]
dL <- reshape2::melt(dW,
                     id.vars = c("G","Id","Gender","X1","X2"),
                     variable.name = "time")
dL <- dL[order(dL$Id),,drop=FALSE]

## *** model fit
m <- lvm(c(Y1[mu1:sigma]~1*eta,
           Y2[mu2:sigma]~1*eta,
           Y3[mu3:sigma]~1*eta,
           eta~G+Gender))
e.lvm <- estimate(m, dW)


e.lmer <- lme4::lmer(value ~ time + G + Gender + (1|Id),
                     data = dL, REML = FALSE)


## *** tests
Sigma0.lvm <- vcov(e.lvm)
SigmaAdj.lvm <- attr(residuals2(e.lvm, bias.correct = TRUE, return.vcov.param = TRUE), "vcov.param")
## Sigma0.lvm/SigmaAdj.lvm
SigmaAdjRed.lvm <- SigmaAdj.lvm[1:5,1:5]
Sigma.GS <- vcovAdj(e.lmer)

## not far from KR correction
SigmaAdjRed.lvm/Sigma.GS
Sigma0.lvm[1:5,1:5]/Sigma.GS

## * leverage adjusted score for lvm models
m.sim <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m.sim) <- eta1~X1+X2
latent(m.sim) <- ~eta1+eta2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)

## ** linear regression
e0 <- estimate(lvm(Y1~X1),d)
s0 <- score2(e0, bias.correct = TRUE)

test_that("score2.lvm vs score2.lm (adj)",{
    sGS <- score2(lm(Y1~X1, data = d), bias.correct = TRUE)
    expect_equal(unname(s0),unname(sGS))
})

## ** multiple linear regression
e1 <- estimate(lvm(Y1~X1,Y2~X2+X3,Y3~1),d)

test_that("score2.lvm(adj): univariate vs. multiple univariate",{
    s1 <- score2(e1, bias.correct = TRUE)
    expect_equal(s1[,c("Y1","Y1~X1","Y1~~Y1")],
                 s0[,c("Y1","Y1~X1","Y1~~Y1")])
})



## ** lvm model with leverage adjusted residuals
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1

e <- estimate(m,d)

test_that("",{
    r2 <- score2(e)
})

##----------------------------------------------------------------------
### test1-sCorrect-smallSampleCorrection.R ends here

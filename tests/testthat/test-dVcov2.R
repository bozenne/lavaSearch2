### test-dVcov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2018 (15:17) 
## Version: 
## Last-Updated: jan  4 2018 (15:05) 
##           By: Brice Ozenne
##     Update #: 23
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

context("dVcov2")
n <- 5e1

mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G+Gender,X1~1,X2~1))
latent(mSim) <- ~eta
categorical(mSim, labels = c("M","F")) <- ~Gender
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- as.data.table(sim(mSim,n,latent = FALSE))
setkey(dW, "Id")
dL <- melt(dW,id.vars = c("G","Id","Gender","X1","X2"), variable.name = "time")
setkey(dL, "Id")

## * linear regression

## ** lm
e.lvm <- estimate(lvm(Y1~X1+X2), data = dW)
e.lvm$prepareScore2 <- prepareScore2(e.lvm, second.order = TRUE, usefit = FALSE)
e.gls <- gls(Y1~X1+X2, data = dW, method = "ML")

test_that("linear regression: dVcov2",{
    ## lvm
    GS.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lvm, res.lvm)

    ## gls
    GS.gls <- dVcov2(e.gls, cluster = dW$Id, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.gls <- dVcov2(e.gls, cluster = dW$Id, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.gls, res.gls)
    
})

## * mixed model

## ** Compound symmetry
m <- lvm(c(Y1[mu1:sigma]~1*eta,
           Y2[mu2:sigma]~1*eta,
           Y3[mu3:sigma]~1*eta,
           eta~G+Gender)) 
e.lvm <- estimate(m, dW)

e.lmer <- lmer(value ~ time + G + Gender + (1|Id),
               data = dL, REML = FALSE)

e.lme <- lme(value ~ time + G + Gender, random = ~ 1|Id, data = dL, method = "ML")
e.gls <- gls(value ~ time + G + Gender,
             correlation = corCompSymm(form=~ 1|Id),
             data = dL, method = "ML")

expect_equal(as.double(logLik(e.lmer)),as.double(logLik(e.lvm)))

test_that("mixed model CS: dVcov2",{
    ## lvm
    GS.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lvm, res.lvm)

    ## gls
    GS.gls <- dVcov2(e.gls, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.gls <- dVcov2(e.gls, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.gls, res.gls)

    ## lme
    GS.lme <- dVcov2(e.lme, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lme <- dVcov2(e.lme, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lme, res.lme)
})

## ** Unstructured with weights

m <- lvm(c(Y1~eta,Y2~eta,Y3~eta,eta~G+Gender))
e.lvm <- estimate(m, dW)
e.lme <- lme(value ~ time + G + Gender,
             random = ~ 1|Id,
             correlation = corSymm(),
             weight = varIdent(form = ~1|time),
             data = dL, method = "ML")
e.gls <- gls(value ~ time + G + Gender,
             correlation = corSymm(form=~ 1|Id),
             weight = varIdent(form = ~1|time),
             data = dL, method = "ML")
e.gls <- gls(value ~ 1,#time + G + Gender,
             correlation = corSymm(form=~ 1|Id),
             weight = varIdent(form = ~1|time),
             data = dL, method = "ML")

logLik(e.lvm)
logLik(e.lme)
logLik(e.gls)

test_that("mixed model UN: dVcov2",{
    ## lvm
    GS.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lvm, res.lvm)

    ## gls
    GS.gls <- dVcov2(e.gls, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.gls <- dVcov2(e.gls, adjust.residuals = FALSE,
                      numericDerivative = FALSE)

    expect_equal(GS.gls, res.gls)

    ## lme
    ## pb: singular information matrix
    ## GS.lme <- dVcov2(e.lme, adjust.residuals = FALSE,
    ##                  numericDerivative = TRUE)
    ## res.lme <- dVcov2(e.lme, adjust.residuals = FALSE,
    ##                   numericDerivative = FALSE)
})

## * latent variable model

m.sim <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
               Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m.sim) <- eta1~X1+X2
latent(m.sim) <- ~eta1+eta2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)


## ** 1 factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
latent(m) <- ~eta1
regression(m) <- eta1~X1+X2

m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1))
latent(m) <- ~eta1

e.lvm1F <- estimate(m,d)

GS.lvm1F <- dVcov2(e.lvm1F, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
test_that("1 factor model (at ML)",{    
    res.lvm1F <- dVcov2(e.lvm1F, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lvm1F, res.lvm1F)
    ##    dfVariance(e.lvm1F)
})

## ** 2 factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1+eta2

m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1,
           Z1~eta2,Z2~eta2,Z3~eta2))
latent(m) <- ~eta1+eta2

e.lvm2F <- estimate(m,d)

GS.lvm2F <- dVcov2(e.lvm2F, adjust.residuals = FALSE,
                   numericDerivative = TRUE)
test_that("2 factor model (at ML)",{
    res.lvm2F <- dVcov2(e.lvm2F, adjust.residuals = FALSE,
                        numericDerivative = FALSE)
    
    expect_equal(GS.lvm2F, res.lvm2F)

    expect_equal(GS.lvm2F[dimnames(GS.lvm1F)[[1]],dimnames(GS.lvm1F)[[2]],dimnames(GS.lvm1F)[[3]]],
             GS.lvm1F)
    ##    dfVariance(e.lvm)
})

expect_equal(res.lvm2F[dimnames(GS.lvm1F)[[1]],dimnames(GS.lvm1F)[[2]],dimnames(GS.lvm1F)[[3]]],
             GS.lvm1F)

round(res.lvm2F[dimnames(GS.lvm1F)[[1]],dimnames(GS.lvm1F)[[2]],"eta1"]-GS.lvm1F[,,"eta1"],7)

attr(GS.lvm1F, "vcov.param") <- NULL
attr(GS.lvm1F, "param") <- NULL
attr(GS.lvm2F, "vcov.param") <- NULL
attr(GS.lvm2F, "param") <- NULL
attr(res.lvm2F, "vcov.param") <- NULL
attr(res.lvm2F, "param") <- NULL

## ** 2 factor model (covariance)

## ** 2 factor model (correlation)

##----------------------------------------------------------------------
### test-dVcov2.R ends here

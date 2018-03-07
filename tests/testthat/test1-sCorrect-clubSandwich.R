### test1-sCorrect-clubSandwich.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  7 2018 (12:08) 
## Version: 
## Last-Updated: mar  7 2018 (18:13) 
##           By: Brice Ozenne
##     Update #: 16
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
library(clubSandwich)
library(nlme)
context("sCorrect: replicate clubSandwich results")

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
dL <- reshape2::melt(d, id.vars = c("Id","X1","X2","X3","Gender"),
                     measure.vars = c("Y1","Y2","Y3","Z1","Z2","Z3"))
dLred <- dL[dL$variable %in% c("Y1","Y2","Y3"),]
dLred$variable.factor <- as.factor(dLred$variable)

## * linear regression [lm,gls,lvm]
## ** model fit and sCorrect
e.lm <- lm(Y1~X1+X2, data = d)

## ** iid2 matches clubSandwich
test_that("iid2.lm/iid2.lvm matches clubSandwich", {
    V.GS <- clubSandwich::vcovCR(e.lm, type = "CR2", cluster = d$Id)

    eHC2.iid2.lm <- iid2(e.lm, value = TRUE)
    V.lm <- crossprod(eHC2.iid2.lm)

    expect_equal(as.matrix(V.GS),
                 V.lm[rownames(V.GS),colnames(V.GS)],
                 tol = 1e-7)
    ## a <- sCorrect(e.lm,adjust.n=TRUE,adjust.Omega=TRUE, n.iter = 1)$Omega 
    ## b <- sCorrect(e.lm,adjust.n=FALSE,adjust.Omega=FALSE)$Omega
    ## a/b
})


## * multiple linear regression with constrains [lvm, gls]
## ** model fit and sCorrect
e.gls <- gls(value ~ X1 + X2,
             data = dLred,
             weight = varIdent(form = ~1|variable),
             method = "ML")
e.lvm <- estimate(lvm(Y1[mu:sigma1]~ beta1*X1 + beta2*X2,
                      Y2[mu:sigma2]~ beta1*X1 + beta2*X2,
                      Y3[mu:sigma3]~ beta1*X1 + beta2*X2),
                  data = d)

factor <- (e.gls$dims$N - e.gls$dims$p)/(e.gls$dims$N - e.gls$dims$p * (e.gls$method == "REML"))
index.coef <- 1:3

test_that("gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.gls)))
})
 
## ** HC0/HC1
iid2HC0.gls <- iid2(e.gls, value = FALSE, cluster = "Id")
iid2HC0.lvm <- iid2(e.lvm, value = FALSE)

test_that("gls: HC0/HC1", {
    expect_equal(unname(iid2HC0.gls[,index.coef]),
                 unname(iid2HC0.lvm[,index.coef]),
                 tol = 1e-5)
    
    VHC0.gls <- crossprod(iid2HC0.gls)[index.coef,index.coef]
    GS <- clubSandwich::vcovCR(e.gls, type = "CR0", cluster = dLred$Id) * factor^2
    expect_equal(as.double(GS),as.double(VHC0.gls), tolerance = 1e-10)
    
    GS <- clubSandwich::vcovCR(e.gls, type = "CR1", cluster = dLred$Id) * factor^2
    VHC1.gls <- VHC0.gls*n/(n-1)
    expect_equal(as.double(GS),as.double(VHC1.gls), tolerance = 1e-10)
})

## ** HC2
iid2HC2.gls <- iid2(e.gls, value = TRUE, n.iter = 1, cluster = "Id")
iid2HC2.lvm <- iid2(e.lvm, value = TRUE, n.iter = 1)

test_that("gls: HC2", {
    expect_equal(unname(iid2HC2.gls[,index.coef]),
                 unname(iid2HC2.lvm[,index.coef]),
                 tol = 1e-5)

    VHC2.gls <- crossprod(iid2HC2.gls)[index.coef,index.coef]
    GS <- clubSandwich::vcovCR(e.gls, type = "CR2", cluster = dLred$Id) * factor^2
    expect_equal(as.double(GS),as.double(VHC2.gls), tolerance = 1e-10)
})

## * mixed model: CS [lvm,gls,lme]
## ** model fit and sCorrect
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
index.coef <- 1:length(coef(e.gls))

test_that("lme/gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lme)))
})

factor <- (e.gls$dims$N - e.gls$dims$p)/(e.gls$dims$N - e.gls$dims$p * (e.gls$method == "REML"))

## ** HC0/HC1
iid2HC0.gls <- iid2(e.gls, value = FALSE)
iid2HC0.lme <- iid2(e.lme, value = FALSE)
iid2HC0.lvm <- iid2(e.lvm, value = FALSE)

test_that("gls: HC0/HC1", {
    expect_equal(unname(iid2HC0.gls[,index.coef]),
                 unname(iid2HC0.lvm[,index.coef]),
                 tol = 1e-5)
    
    VHC0.gls <- crossprod(iid2HC0.gls)[index.coef,index.coef]
    GS <- clubSandwich::vcovCR(e.gls, type = "CR0", cluster = dLred$Id) * factor^2
    expect_equal(as.double(GS),as.double(VHC0.gls), tolerance = 1e-10)

    GS <- clubSandwich::vcovCR(e.gls, type = "CR1", cluster = dLred$Id) * factor^2
    VHC1.gls <- VHC0.gls*n/(n-1)
    expect_equal(as.double(GS),as.double(VHC1.gls), tolerance = 1e-10)
})

## ** HC2
iid2HC2.gls <- iid2(e.gls, value = TRUE, n.iter = 1)
iid2HC2.lme <- iid2(e.lme, value = TRUE, n.iter = 1)
iid2HC2.lvm <- iid2(e.lvm, value = TRUE, n.iter = 1)

test_that("gls: HC2", {
    expect_equal(unname(iid2HC2.gls[,index.coef]),
                 unname(iid2HC2.lvm[,index.coef]),
                 tol = 1e-5)

    VHC2.gls <- crossprod(iid2HC2.gls)[index.coef,index.coef]
    GS <- clubSandwich::vcovCR(e.gls, type = "CR2", cluster = dLred$Id) * factor^2
    expect_equal(as.double(GS),as.double(VHC2.gls), tolerance = 1e-10)
})


## * mixed model: CS with different variances [lvm,gls,lme]
m <- lvm(c(Y1[0:sigma1]~1*eta,
           Y2[0:sigma2]~1*eta,
           Y3[0:sigma3]~1*eta,
           eta~X1+X2))
latent(m) <- ~eta
e.lvm <- estimate(m, d)

e.lme <- nlme::lme(value ~ X1 + X2,
                   random =~1| Id,
                   weights = varIdent(form =~ 1|variable),
                   data = dLred, method = "ML")

e.gls <- nlme::gls(value ~ X1 + X2,
                   correlation = corCompSymm(form = ~1| Id),
                   weights = varIdent(form =~ 1|variable),
                   data = dLred, method = "ML")
index.coef <- 1:length(coef(e.gls))

test_that("lme/gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
    ## gls does not give the same likelihood
    ## expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lme)))
})

factor <- (e.gls$dims$N - e.gls$dims$p)/(e.gls$dims$N - e.gls$dims$p * (e.gls$method == "REML"))

## ** HC0/HC1
iid2HC0.gls <- iid2(e.gls, value = FALSE)
iid2HC0.lme <- iid2(e.lme, value = FALSE)
iid2HC0.lvm <- iid2(e.lvm, value = FALSE)

iid2HC0.lvm-iid2HC0.lme
coef(e.lvm)[1:3]-fixef(e.lme)

e2.lme <- e.lme
sCorrect(e2.lme) <- FALSE

e2.gls <- e.gls
sCorrect(e2.gls) <- FALSE

e2.lvm <- e.lvm
sCorrect(e2.lvm) <- FALSE

e2.lvm$sCorrect$Omega
e2.lme$sCorrect$Omega
coef()

iid2HC0.lvm[1:5,1:3]-iid2HC0.lme[1:5,1:3]
## BUG IN sCorrect.lme ???

test_that("gls: HC0/HC1", {
    expect_equal(unname(iid2HC0.lme[,index.coef]),
                 unname(iid2HC0.lvm[,index.coef]),
                 tol = 1e-5)

        VHC0.lvm <- crossprod(iid2HC0.lvm)[1:3,1:3]

    VHC0.lme <- crossprod(iid2HC0.lme)[index.coef,index.coef]
    GS <- clubSandwich::vcovCR(e.lme, type = "CR0", cluster = dLred$Id) * factor^2
    expect_equal(as.double(GS),as.double(VHC0.lme), tolerance = 1e-10)

    GS <- clubSandwich::vcovCR(e.gls, type = "CR1", cluster = dLred$Id) * factor^2
    VHC1.gls <- VHC0.gls*n/(n-1)
    expect_equal(as.double(GS),as.double(VHC1.gls), tolerance = 1e-10)
})

## ** HC2
iid2HC2.gls <- iid2(e.gls, value = TRUE, n.iter = 1)
iid2HC2.lme <- iid2(e.lme, value = TRUE, n.iter = 1)
iid2HC2.lvm <- iid2(e.lvm, value = TRUE, n.iter = 1)

test_that("gls: HC2", {
    expect_equal(unname(iid2HC2.gls[,index.coef]),
                 unname(iid2HC2.lvm[,index.coef]),
                 tol = 1e-5)

    VHC2.gls <- crossprod(iid2HC2.gls)[index.coef,index.coef]
    GS <- clubSandwich::vcovCR(e.gls, type = "CR2", cluster = dLred$Id) * factor^2
    expect_equal(as.double(GS),as.double(VHC2.gls), tolerance = 1e-10)
})


test_that("lme/gls/lvm: HC0/HC1", {
    iid2HC0.lme <- iid2(e.lme, data = dL, bias.correct = FALSE)
    iid2HC0.lvm <- iid2(e.lvm, bias.correct = FALSE)    
    expect_equal(unname(iid2HC0.lme[,index.coef]),unname(iid2HC0.lvm[,index.coef]), tol = 1e-6)

    VsandwichHC0.lme <- crossprod(iid2HC0.lme)[index.coef,index.coef]
    VsandwichHC0.lvm <- crossprod(iid2HC0.lvm)[index.coef,index.coef]
    expect_equal(as.double(VsandwichHC0.lme), as.double(VsandwichHC0.lvm), tol = 1e-7)
    GS <- clubSandwich::vcovCR(e.lme, type = "CR0", cluster = dL$Id)
    expect_equal(as.double(GS),as.double(VsandwichHC0.lme))
    
    GS <- clubSandwich::vcovCR(e.lme, type = "CR1", cluster = dL$Id)
    VsandwichHC1.lme <-  VsandwichHC0.lme*n/(n-1)
    expect_equal(as.double(GS),as.double(VsandwichHC1.lme))    
})

test_that("lme/gls/lvm: HC2", {
    iid2HC2.lme <- iid2(e.lme, data = dL, bias.correct = TRUE, as.clubSandwich = 2)
    iid2HC2.lvm <- iid2(e.lvm, bias.correct = TRUE, as.clubSandwich = 2)    
    expect_equal(unname(iid2HC2.lme[,index.coef]),unname(iid2HC2.lvm[,index.coef]), tol = 1e-6)

    VsandwichHC2.lme <- crossprod(iid2HC2.lme)[index.coef,index.coef]
    VsandwichHC2.lvm <- crossprod(iid2HC2.lvm)[index.coef,index.coef]
    expect_equal(as.double(VsandwichHC2.lme), as.double(VsandwichHC2.lvm), tol = 1e-7)
    GS <- clubSandwich::vcovCR(e.lme, type = "CR2", cluster = dL$Id)
    expect_equal(as.double(GS),as.double(VsandwichHC2.lme))
})

## * gls - Unstructured
m <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~1))
covariance(m) <- Y1~Y2
covariance(m) <- Y1~Y3
e.lvm <- estimate(m, dW)

dataRed <- dL[dL$time!="Y4",]
dataRed$time <- droplevels(dataRed$time)
e.gls <- nlme::gls(value ~ time,
                   correlation = corSymm(form =~ 1| Id),
                   weight = varIdent(form =~ 1|time),
                   data = dataRed, method = "ML")

e.lme <- nlme::lme(value ~ time,
                   random =~ 1|Id,
                   correlation = corSymm(),
                   weight = varIdent(form =~ 1|time),
                   data = dataRed, method = "ML")

index.coef <- 1:length(coef(e.gls))

test_that("lme/gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lme)))
})


test_that("gls/lvm: HC0/HC1", {
    iid2HC0.gls <- iid2(e.gls, data = dataRed, bias.correct = FALSE)
    ## iid2HC0.lme <- iid2(e.lme, bias.correct = FALSE) ## not invertible
    iid2HC0.lvm <- iid2(e.lvm, bias.correct = FALSE)
    

    VsandwichHC0.gls <- crossprod(iid2HC0.gls)[index.coef,index.coef]
    VsandwichHC0.lvm <- crossprod(iid2HC0.lvm)[index.coef,index.coef]    
    expect_equal(as.double(VsandwichHC0.gls), as.double(VsandwichHC0.lvm), tol = 1e-7)
    ## scoreX <- score2(e.lvm, bias.correct = FALSE)
    ## colSums(scoreX)
    GS <- clubSandwich::vcovCR(e.gls, type = "CR0", cluster = dataRed$Id)
    ## VsandwichHC0.gls-GS
    ## VsandwichHC0.lvm-GS

    ## VsandwichHC0.gls-vcov(e.gls)
    ## GS-vcov(e.gls)

})

test_that("lme: HC2", {
    iid2HC2.gls <- iid2(e.gls, data = dataRed, bias.correct = TRUE, as.clubSandwich = 2)
    iid2HC2.lvm <- iid2(e.lvm, bias.correct = TRUE, as.clubSandwich = 2)    

    VsandwichHC2.gls <- crossprod(iid2HC2.gls)[index.coef,index.coef]
    VsandwichHC2.lvm <- crossprod(iid2HC2.lvm)[index.coef,index.coef]    
    expect_equal(as.double(VsandwichHC2.gls), as.double(VsandwichHC2.lvm), tol = 1e-7)

    ## VsandwichHC2.gls <- crossprod(iid2HC2.gls)[index.coef,index.coef]
    ## GS <- clubSandwich::vcovCR(e.gls, type = "CR2", cluster = dL$Id)
    ## expect_equal(as.double(GS),as.double(VsandwichHC2.gls))
})


##----------------------------------------------------------------------
### test1-sCorrect-clubSandwich.R ends here

### test1-sCorrect-compatibilityWithLava.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  6 2018 (10:40) 
## Version: 
## Last-Updated: mar  6 2018 (14:00) 
##           By: Brice Ozenne
##     Update #: 21
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## source("c:/Users/hpl802/Documents/GitHub/lavaSearch2/tests/testthat/test1-iid2-lava.R")

## * header
if(FALSE){ ## already called in test-all.R
    rm(list = ls())
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))
library(nlme)
context("sCorrect: replicate lava results")

## * linear regression
n <- 5e1
m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e.lvm <- estimate(m,d)
e.lm <- lm(Y~X1+X2+X3, data = d)
e.gls <- gls(Y~X1+X2+X3, data = d, method = "ML")
param <- coef(e.lvm)
e2.lvm <- e.lvm
e2.gls <- e.gls
e2.lm <- e.lm

## ** sCorrect
sCorrect(e2.lvm) <- FALSE
sCorrect(e2.gls, cluster = 1:n) <- FALSE
sCorrect(e2.lm) <- FALSE

## ** check score, iid, residuals, vcov, compare2 at ML
test_that("linear regression (at ML) internal consistency",{
    expect_equivalent(e2.lvm$sCorrect$Omega,e2.lm$sCorrect$Omega)
    expect_equivalent(e2.gls$sCorrect$Omega,e2.lm$sCorrect$Omega)

    expect_equivalent(e2.lvm$sCorrect$vcov.param,e2.lm$sCorrect$vcov.param)
    expect_equivalent(e2.gls$sCorrect$vcov.param,e2.lm$sCorrect$vcov.param)

    expect_equivalent(e2.lvm$sCorrect$leverage,e2.lm$sCorrect$leverage)
    expect_equivalent(e2.gls$sCorrect$leverage,e2.lm$sCorrect$leverage)

    expect_equivalent(e2.lvm$sCorrect$score,e2.lm$sCorrect$score)
    expect_equivalent(e2.gls$sCorrect$score,e2.lm$sCorrect$score)

    expect_equivalent(e2.lvm$sCorrect$epsilon,e2.lm$sCorrect$epsilon)
    expect_equivalent(e2.gls$sCorrect$epsilon,e2.lm$sCorrect$epsilon)

    expect_equivalent(e2.lvm$sCorrect$dVcov.param,e2.lm$sCorrect$dVcov.param)
    expect_equivalent(e2.gls$sCorrect$dVcov.param,e2.lm$sCorrect$dVcov.param)    
}

test_that("linear regression (at ML) compare to lava",{

    expect_equivalent(e2.lvm$sCorrect$vcov.param, vcov(e.lvm))
    expect_true(all(e2.lvm$sCorrect$leverage==0))
    expect_true(e2.lvm$sCorrect$n.corrected==e.lvm$data$n)
    expect_equal(e2.lvm$sCorrect$epsilon, residuals(e.lvm))
    expect_equal(e2.lvm$sCorrect$param, coef(e.lvm))
    expect_equal(e2.lvm$sCorrect$score, score(e.lvm, indiv = TRUE))

    GS <- iid(e.lm)
    expect_equivalent(iid2(e2.lvm)[,1:length(colnames(GS))], GS)
    expect_equal(e2.lvm$sCorrect$epsilon, residuals2(e2.lvm))
    expect_true(all(leverage2(e2.lvm) == 0))

    C1 <- compare2(e2.lvm, par = c("Y~X1","Y~X2"))
    C2 <- lava::compare(e.lvm, par = c("Y~X1","Y~X2"))
    expect_equal(unname(C1$statistic),
                 unname(C2$statistic/NROW(C1$estimate))
                 )
                 
    
})

## ** check score not at ML
test_that("linear regression (at ML + 1) compare to lava",{
    S1 <- score2(e2.lvm, param = coef(e2.lvm)+1)
    S2 <- score2(e.lvm, param = coef(e.lvm)+1, value = FALSE)

    S3 <- score2(e2.lm, param = lavaSearch2:::.coef2(e2.gls)+1) ## not .coef2(e2.lm) because different estimate of the variance
    S4 <- score2(e.lm, param = lavaSearch2:::.coef2(e.gls)+1, value = FALSE)

    S5 <- score2(e2.gls, param = lavaSearch2:::.coef2(e2.gls)+1)
    S6 <- score2(e.gls, param = lavaSearch2:::.coef2(e.gls)+1, cluster = 1:n, value = FALSE)

    GS <- score(e.lvm, p = coef(e.lvm)+1, indiv = TRUE)

    
    expect_equal(S1, GS)
    expect_equal(S2, GS)
    
    expect_equal(as.double(S3), as.double(GS))
    expect_equal(as.double(S4), as.double(GS))
    expect_equal(as.double(S5), as.double(GS))
    expect_equal(as.double(S6), as.double(GS))
    
})


test_that("linear regression (at ML + 1:p) compare to lava",{
    p <- length(coef(e.lvm))
    S1 <- score2(e2.lvm, param = coef(e2.lvm)+1:p)
    S2 <- score2(e.lvm, param = coef(e.lvm)+1:p, value = FALSE)

    S3 <- score2(e2.lm, param = lavaSearch2:::.coef2(e2.gls)+1:p) ## not .coef2(e2.lm) because different estimate of the variance
    S4 <- score2(e.lm, param = lavaSearch2:::.coef2(e.gls)+1:p, value = FALSE)

    S5 <- score2(e2.gls, param = lavaSearch2:::.coef2(e2.gls)+1:p)
    S6 <- score2(e.gls, param = lavaSearch2:::.coef2(e.gls)+1:p, cluster = 1:n, value = FALSE)

    GS <- score(e.lvm, p = coef(e.lvm)+1:p, indiv = TRUE)
    
    expect_equal(S1, GS)
    expect_equal(S2, GS)
    
    expect_equal(as.double(S3), as.double(GS))
    expect_equal(as.double(S4), as.double(GS))
    expect_equal(as.double(S5), as.double(GS))
    expect_equal(as.double(S6), as.double(GS))
    
})


## ** check score with constrains

test_that("linear regression: constrains",{
    m <- lvm(Y[0:2]~X1+1*X2)
    e <- estimate(m, d)
    expect_equal(score2(e, bias.correct = FALSE),
                 score(e, indiv = TRUE))
    
    m <- lvm(Y~beta*X1+beta*X2)
    e <- estimate(m, d)
    expect_equal(score2(e, bias.correct = FALSE),
                 score(e, indiv = TRUE))
})


## * multiple linear regression
## ** without covariance link
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)
## prepareScore2(e) <- FALSE

test_that("multiple linear regression (at ML)",{
    test <- score2(e, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("multiple linear regression (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("multiple linear regression (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

})

test_that("multiple linear regressions: constrains",{
    m <- lvm(Y1~X1+1*X2,Y2~2*X3+2*X1,Y3~X2)
    e <- estimate(m, d)

    expect_equal(score2(e, bias.correct = FALSE),
                 score(e, indiv = TRUE))
})

## ** with covariance links
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
covariance(m) <- Y1~Y2
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)
## prepareScore2(e) <- FALSE

test_that("multiple linear regression, covariance link (at ML)",{
    test <- score2(e, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, indiv=TRUE)    
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("multiple linear regression, covariance link (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("multiple linear regression, covariance link (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)
})

## * latent variable model
m.sim <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m.sim) <- eta1~X1+X2
latent(m.sim) <- ~eta1+eta2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)

## ** factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
regression(m) <- eta1~X1+X2

e <- estimate(m,d)
param <- coef(e)
## prepareScore2(e) <- FALSE

test_that("factor model (at ML)",{
    test <- score2(e, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
     
    test <- score2(e, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})
test_that("factor model (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("factor model (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

})

test_that("factor model: fixed coefficients",{
    m <- lvm(Y1~1*eta+1*X2,Y2~1*eta,Y3~1*eta)
    e <- estimate(m, d)

    expect_equal(score2(e, bias.correct = FALSE),
                 score(e, indiv = TRUE))
})

test_that("factor model: constrains",{
    m <- lvm(Y1~1*eta+X2,Y2~lambda*eta+X2,Y3~lambda*eta,eta ~ beta*X2+beta*X1)
    e <- estimate(m, d)

    expect_equal(score2(e, bias.correct = FALSE),
                 score(e, indiv = TRUE))
})


## ** 2 factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
param <- coef(e)
## prepareScore2(e) <- FALSE

test_that("2 factor model (at ML)",{
    test <- score2(e, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

})

test_that("2 factor model: constrains",{
    m <- lvm(Y1~1*eta1+X2,Y2~lambda*eta1+X2,Y3~lambda*eta1,eta1 ~ beta*X2+beta*X1,
             Z1~0+eta2,Z2~lambda*eta2,Z3~eta2)
    e <- estimate(m, d)

    expect_equal(score2(e, bias.correct = FALSE),
                 score(e, indiv = TRUE))
})

## ** 2 factor model (covariance)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
covariance(m) <- eta1 ~ eta2
## covariance(m) <- Y1 ~ Z1
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
param <- coef(e)
## prepareScore2(e) <- FALSE

test_that("2 factor model, covariance (at ML)",{
    test <- score2(e, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model, covariance (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model, covariance (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

})

## ** 2 factor model (correlation LV)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta2 ~ X1
regression(m) <- eta1 ~ eta2+X2+X3

e <- estimate(m,d)
param <- coef(e)
## prepareScore2(e) <- FALSE

test_that("2 factor model, correlation LV (at ML)",{
    test <- score2(e, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model, correlation LV (not at ML: +1)",{
    test <- score2(e, p = param+1, indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, bias.correct = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("2 factor model, covariancecorrelation LV (not at ML: +1:p)",{
    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, bias.correct = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

})


##----------------------------------------------------------------------
### test1-sCorrect-compatibilityWithLava.R ends here

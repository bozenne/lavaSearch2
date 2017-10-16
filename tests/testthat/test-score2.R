### test-score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 13 2017 (11:28) 
## Version: 
## last-updated: okt 16 2017 (18:33) 
##           By: Brice Ozenne
##     Update #: 30
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(testthat)

context("score2")
n <- 5e1

## * linear regression

m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
test_that("",{
    test <- score2(e)
    GS <- score(e,indiv=TRUE)
    expect_equal(dim(test),dim(GS))
    expect_equal(test[,colnames(GS)],GS)
})

## * multiple linear regression

## ** without covariance link
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
test_that("",{
    test <- score2(e)
    GS <- score(e,indiv=TRUE)
    expect_equal(dim(test),dim(GS))
    expect_equal(test[,colnames(GS)],GS)
})

## ** without covariance links
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
covariance(m) <- Y1~Y2
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
test_that("",{
    test <- score2(e)
    GS <- score(e,indiv=TRUE)
    expect_equal(dim(test),dim(GS))
    expect_equal(test[,colnames(GS)],GS)
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
latent(m) <- ~eta1

e <- estimate(m,d)
test_that("",{
    test <- score2(e)
    GS <- score(e,indiv=TRUE)
    expect_equal(dim(test),dim(GS))
    expect_equal(test[,colnames(GS)],GS)
})

## ** 2 factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
test_that("",{
    test <- score2(e)
    GS <- score(e,indiv=TRUE)
    expect_equal(dim(test),dim(GS))
    expect_equal(test[,colnames(GS)],GS)
})

## ** 2 factor model (covariance)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
covariance(m) <- eta1 ~ eta2
latent(m) <- ~eta1+eta2
e <- estimate(m,d)

test_that("",{
    test <- score2(e)
    GS <- score(e,indiv=TRUE)
    expect_equal(dim(test),dim(GS))
    expect_equal(test[,colnames(GS)],GS)
})

## ** 2 factor model (correlation LV)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta2 ~ X1
regression(m) <- eta1 ~ eta2+X2

## m <- lvm(c(Y1~0+1*eta1,Y2~0+1*eta1,Y3~0+1*eta1,
##            Z1~0+1*eta2,Z2~0+1*eta2,Z3~0+1*eta2))
## regression(m) <- eta2 ~ 0
## regression(m) <- eta1 ~ 0+eta2
## latent(m) <- ~ eta1 + eta2

e <- estimate(m,d)

test_that("",{
    test <- score2(e)
    GS <- score(e,indiv=TRUE)
    expect_equal(dim(test),dim(GS))
    expect_equal(test[,colnames(GS)],GS)
})

## > head(GS)[,"eta1~eta2"]
## [1]  0.78165615 -0.25528969 -1.77126130 -0.36471102 -0.04859673  0.93779653

#----------------------------------------------------------------------
### test-score2.R ends here

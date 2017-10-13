### test-score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 13 2017 (11:28) 
## Version: 
## last-updated: okt 13 2017 (16:04) 
##           By: Brice Ozenne
##     Update #: 14
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
## > GS[,"eta1~~eta1"]
##  [1] -0.286369552 -0.449185901 -0.065588591  1.664995240 -0.452587762 -0.271317267 -0.194578318 -0.448190671  0.073127901 -0.442647841 -0.005928534 -0.386910935 -0.011210398 -0.425451802
## [15]  0.732905115 -0.065541944  0.157563373 -0.320305676  0.508661083  0.039943793 -0.444728604  1.485893082  0.034360262  0.183589276  0.760281243  0.036079345  0.614895186 -0.197295750
## [29] -0.208387568  0.166403845 -0.303540677 -0.449267189 -0.453125141 -0.275213831 -0.451595623 -0.437485206 -0.153324317  1.162168045 -0.178238407 -0.445457110 -0.403239726 -0.352371190
## [43] -0.062117361  0.972539214 -0.436848109 -0.121321244  0.690613533 -0.245415238 -0.411971421  0.572761385
##           Y1        Y2        Y3
## Y1 1.0000000 1.0414810 0.8661981
## Y2 1.0414810 1.0846826 0.9021288
## Y3 0.8661981 0.9021288 0.7502992
## [1] -0.4531767
##            [,1]        [,2]        [,3]
## 1   0.548657600 -0.83537278 -0.29580569
## 2  -0.048802734  0.26177352 -0.36154478
## 3  -0.123955981 -0.28255542  1.49927943

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

## ** 2 factor model (complex)
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

head(test[,colnames(GS)]-GS)
round(colMeans(abs(test[,colnames(GS)]-GS)),4)

#----------------------------------------------------------------------
### test-score2.R ends here

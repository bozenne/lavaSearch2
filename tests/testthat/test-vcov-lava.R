### test-vcov.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: nov  6 2017 (11:44) 
## Version: 
## last-updated: nov 20 2017 (16:29) 
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

## * Model vcov
## ** function
rmAttr <- function(x, name.rm = NULL, name.keep){
    if(is.null(name.rm)){
        name.rm <- names(attributes(x))
    }
    for(iAttr in name.rm){
        attr(x, iAttr) <- NULL
    }
    return(x)
}

## ** linear regression
m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- d

test_that("linear regression (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("linear regression (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("linear regression (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

test_that("linear regression: constrains",{
    m <- lvm(Y[0:2]~X1+1*X2)
    e <- estimate(m, d)

    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))  
    
    m <- lvm(Y~beta*X1+beta*X2)
    e <- estimate(m, d)

    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))  
})


## ** multiple linear regression
## *** without covariance link
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)

prepareScore2(e) <- d

test_that("multiple linear regression (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))

    expect_equal(unname(test), unname(GS))    
})

test_that("multiple linear regression (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("multiple linear regression (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

test_that("multiple linear regressions: constrains",{
    m <- lvm(Y1~X1+1*X2,Y2~2*X3+2*X1,Y3~X2)
    e <- estimate(m, d)

    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

## *** with covariance links
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
covariance(m) <- Y1~Y2
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)
e$prepareScore2 <- prepareScore2(e)

test_that("multiple linear regression, covariance link (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("multiple linear regression, covariance link (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("multiple linear regression, covariance link (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

## ** latent variable model
m.sim <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m.sim) <- eta1~X1+X2
latent(m.sim) <- ~eta1+eta2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)

## *** factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
regression(m) <- eta1~X1+X2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- d

test_that("factor model (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("factor model (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("factor model (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

test_that("factor model: fixed coefficients",{
    m <- lvm(Y1~1*eta+1*X2,Y2~1*eta,Y3~1*eta)
    e <- estimate(m, d)

    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("factor model: constrains",{
    m <- lvm(Y1~1*eta+X2,Y2~lambda*eta+X2,Y3~lambda*eta,eta ~ beta*X2+beta*X1)
    e <- estimate(m, d)
    
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})


## *** 2 factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- d

test_that("2 factor model (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("2 factor model (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("2 factor model (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

test_that("2 factor model: constrains",{
    m <- lvm(Y1~1*eta1+X2,Y2~lambda*eta1+X2,Y3~lambda*eta1,eta1 ~ beta*X2+beta*X1,
             Z1~0+eta2,Z2~lambda*eta2,Z3~eta2)
    e <- estimate(m, d)

     test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

## *** 2 factor model (covariance)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
covariance(m) <- eta1 ~ eta2
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- d

test_that("2 factor model, covariance (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("2 factor model, covariance (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("2 factor model, covariance (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})

## *** 2 factor model (correlation LV)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta2 ~ X1
regression(m) <- eta1 ~ eta2+X2

e <- estimate(m,d)
param <- coef(e)
prepareScore2(e) <- d

test_that("2 factor model, correlation (at ML)",{
    test <- attr(score2(e, p = pars(e), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})

test_that("2 factor model, correlation (not at ML: +1)",{
    test <- attr(score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+1))
    expect_equal(unname(test), GS)
})

test_that("2 factor model, correlation (not at ML: +1:p)",{
    test <- attr(score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.vcov.param = TRUE),
                 "vcov.param")
    GS <- solve(information(e, p = param+0.1*(1:length(param))))
    expect_equal(unname(test), GS)    
})


#----------------------------------------------------------------------
### test-vcov.R ends here


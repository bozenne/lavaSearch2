### test-score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 13 2017 (11:28) 
## Version: 
## last-updated: okt 19 2017 (16:43) 
##           By: Brice Ozenne
##     Update #: 72
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

## * not-adjusted score
## ** linear regression

m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)

test_that("linear regression",{
    ## at ML
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

## ** multiple linear regression
## *** without covariance link
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)

test_that("multiple linear regression",{
    ## at ML
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

## *** without covariance links
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
covariance(m) <- Y1~Y2
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)

test_that("multiple linear regression (covariance link)",{
        ## at ML
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
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
latent(m) <- ~eta1

e <- estimate(m,d)
param <- coef(e)

test_that("factor model",{
        ## at ML
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

## *** 2 factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
param <- coef(e)

test_that("2 factor model",{
        ## at ML
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

## *** 2 factor model (covariance)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
covariance(m) <- eta1 ~ eta2
latent(m) <- ~eta1+eta2

e <- estimate(m,d)
param <- coef(e)

test_that("2 factor model (covariance)",{
        ## at ML
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

## *** 2 factor model (correlation LV)
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
param <- coef(e)

test_that("2 factor model (correlation LV)",{
    ## at ML
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)

    head(round(test-GS,10))
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE, return.df = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

})









## * leverage adjusted score
m.sim <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(m.sim) <- eta1~X1+X2
latent(m.sim) <- ~eta1+eta2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)

## ** linear regression
e0 <- estimate(lvm(Y1~X1),d)
s0 <- score2(e0, adjust.residuals = TRUE)

e1 <- estimate(lvm(Y1~X1,Y2~X2+X3,Y3~1),d)
s1 <- score2(e1, adjust.residuals = TRUE)

test_that("",{
    expect_equal(s1[,c("Y1","Y1~X1","Y1~~Y1")],
                 s0[,c("Y1","Y1~X1","Y1~~Y1")])
    expect_equal(df.residual(e0),attr(s0,"df"))
    ## expect_equal(df.residual(e1),attr(s1,"df"))

    ## df.residual(e1)
    
##    expect_equal(attr(s1, "df"), NROW(d))
})



## ** lvm model with leverage adjusted residuals

m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
regression(m) <- eta1~X1+X2
latent(m) <- ~eta1

e <- estimate(m,d)


e <- estimate(m,d)
param <- coef(e)

test_that("",{
    r2 <- score2(e)
})

#----------------------------------------------------------------------
### test-score2.R ends here



## * error for tobit and multigroup lvm

m.sim <- lvm(Y~X1+X2,G~1)
categorical(m.sim,K=2,label=c("a","b")) <- ~G+X2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)
e <- estimate(list(m.sim,m.sim),data = split(d,d$G))

expect_error(score2(e))


library(lava.tobit)
m.sim <- lvm(Y~X1)
categorical(m.sim,K=2,labels = c("a","b")) <- ~Y
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)
e <- estimate(lvm(Y~X1),data = d)

expect_error(score2(e))

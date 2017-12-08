### test-vcov.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: nov  6 2017 (11:44) 
## Version: 
## last-updated: dec  7 2017 (17:46) 
##           By: Brice Ozenne
##     Update #: 20
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

## * Corrected vcov
n <- 1e4

m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e.lvm <- estimate(lvm(Y~X1+X2+X3),d)
e.lm <- lm(Y~X1+X2+X3, data = d)
Sigma.lvm <- attr(residuals2(e.lvm, return.vcov.param = TRUE), "vcov.param")
Sigma.lm <- attr(residuals2(e.lm, return.vcov.param = TRUE), "vcov.param")

mean.coef <- names(coef(e.lvm))[1:4]
expect_equal(Sigma.test[mean.coef,mean.coef],
             vcov(e.lvm)[mean.coef,mean.coef] * (n+length(mean.coef))/n)

expect_equal(Sigma.test[mean.coef,mean.coef],
             vcov(e.lvm)[mean.coef,mean.coef] * (n+length(mean.coef))/n)






tr(X %*% solve(t(X) %*% X) %*% t(X))

iObs <- 1
lsH <- lapply(1:n, function(iObs){
    X[iObs,,drop=FALSE] %*% solve(t(X) %*% X) %*% t(X[iObs,,drop=FALSE])
})

t(X) %*% X

dVcov.dtheta

2*((1+n.coef/n)*coef(e)["Y~~Y"])^2/n
iXX*(1+n.coef/n)

vcov.param[keep.param,keep.param,drop=FALSE]

diag(solve(t(X) %*% X ))^2 - diag(solve(t(X) %*% X %*% t(X) %*% X))

H <- X %*% solve(t(X) %*% X) %*% t(X)

epsilon <- residuals(e)
sigma_corrected <- mean(sapply(1:n, function(iObs){    
    epsilon[iObs]^2/(1-H[iObs,iObs])
}))
n/(n-4)*coef(e)["Y~~Y"]
coef(e)["Y~~Y"]*(1+4/n)
sigma_corrected


vcov.lm <- attr(residuals2(e.lm, adjust.residuals = TRUE, return.vcov.param = TRUE), "vcov.param")
vcov.lvm <- attr(residuals2(e.lvm, adjust.residuals = TRUE, return.vcov.param = TRUE), "vcov.param")

vcov.lm/vcov(e.lm)
vcov.lvm[1:4,1:4]/vcov(e.lm)

library(pbkrtest)

data(beets, package='pbkrtest')
lg <- lmer(sugpct ~ block + sow + harvest + (1|block:harvest),
data=beets, REML=FALSE)
xx <- KRmodcomp(lg, sm)

getKR(xx)


(fmLarge <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
## removing Days
(fmSmall <- lmer(Reaction ~ 1 + (Days|Subject), sleepstudy))
anova(fmLarge,fmSmall)
KRmodcomp(fmLarge, fmSmall)  ## 17 denominator df




fm1 <- lmer(value ~ G + (1|Id), data = dL)
logLik(fm1)
fm2 <- lme(value ~ G,
           random = ~ 1|Id,
           data = dL)
logLik(fm2)

vcov(fm1)
v1 <- vcovAdj(fm1, detail = 1)
v1/vcov(fm1)

v2 <- attr(residuals2(fm2, adjust.residuals = TRUE, return.vcov.param = TRUE),
           "vcov.param")
v2[1:2,1:2]/vcov(fm1)
## Here the adjusted and unadjusted covariance matrices are identical,
## but that is not generally the case
v1 <- vcov(fm1)
v2 <- vcovAdj(fm1,detail=0)

#----------------------------------------------------------------------
### test-vcov.R ends here





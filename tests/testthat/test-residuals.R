### test-residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (09:08) 
## Version: 
## Last-Updated: nov  9 2017 (13:41) 
##           By: Brice Ozenne
##     Update #: 24
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

context("residuals2")

n <- 5e1

mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
latent(mSim) <- ~eta
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- as.data.table(sim(mSim,n,latent = FALSE))
setkey(dW, "Id")
dL <- melt(dW,id.vars = c("G","Id"), variable.name = "time")
setkey(dL, "Id")
dL$Z1 <- rnorm(NROW(dL))

## * univariate linear model
m <- lvm(Y~X)
d <- sim(m,1e2)

test_that("residuals2 match residuals.lm", {
    GS <- residuals(lm(Y~X, data = d))
    e <- estimate(lvm(Y~X), d)
    res <- residuals(e)
    res2 <- residuals2(e, adjust.residuals = FALSE)

    expect_equal(as.double(res),as.double(GS))
    expect_equal(as.double(res2),as.double(GS))
})


## * multivariate linear models
m <- lvm(Y~G+X,G~X)
d <- sim(m,1e2)

test_that("residuals2 match residuals.lm", {
    GS <- cbind(residuals(lm(Y~1, data = d)),
                residuals(lm(G~1, data = d)))
    e <- estimate(lvm(Y~1,G~1), d)
    res2 <- residuals2(e, adjust.residuals = FALSE)
    res <- residuals(e)

    expect_equal(unname(res),unname(GS))
    expect_equal(unname(res2),unname(GS))

    GS <- cbind(residuals(lm(Y~G+X, data = d)),
                residuals(lm(G~1, data = d)))
    e <- estimate(lvm(Y~G+X,G~1), d)

    res2 <- residuals2(e, adjust.residuals = FALSE)
    res <- residuals(e)

    
    ## expect_equal(as.double(res),as.double(GS))
    ## note: vcov(lm(Y~G+X, data = d))/vcov(e)[c("Y","Y~G","Y~X"),c("Y","Y~G","Y~X")]
    
    expect_equal(unname(coef(e)[c("Y","Y~G","Y~X")]),
                 unname(coef(lm(Y~G+X,data=d))), tol = 1e-5)
    expect_equal(res2[,"Y"],unname(GS[,1]), tol = 1e-4)
    expect_equal(res2[,"G"],unname(GS[,2]))
    
    GS <- cbind(residuals(lm(Y~G+X, data = d)),
                residuals(lm(G~X, data = d)))
    e <- estimate(lvm(Y~G+X,G~X), d)

    res2 <- residuals2(e, adjust.residuals = FALSE)
    res <- residuals(e)

    ## expect_equal(as.double(res),as.double(GS))
    ## note: vcov(lm(Y~G+X, data = d))/vcov(e)[c("Y","Y~G","Y~X"),c("Y","Y~G","Y~X")]
    expect_equal(as.double(res2),as.double(GS))
})

## * mixed model
## ** versus nlme
mSim <- lvm(c(Y1~1*eta1,Y2~1*eta1,Y3~1*eta1,eta1~G1))
latent(mSim) <- ~eta1
transform(mSim, Id~Y1) <- function(x){1:NROW(x)}
dW <- as.data.table(sim(mSim, 5e1, latent = FALSE))
dL <- melt(dW, id.vars = c("Id","G1")) 


test_that("equivalence residuals2.lvm residuals.lvm", {
    m <- lvm(c(Y1[mu1:sigma]~1*eta1,Y2[mu2:sigma]~1*eta1,Y3[mu3:sigma]~1*eta1,eta1~G1))
    latent(m) <- ~eta1
    e.lvm <- estimate(m,dW)

    e.gls <- gls(value ~ variable + G1, data = dL,
                 correlation = corCompSymm(form =~ variable|Id))
    e.lme <- lme(value ~ variable + G1, data = dL,
                 random =~ 1|Id)

    ##   logLik(e.lvm)
    ##   logLik(e.gls)    # not equal!
    ##   logLik(e.lme)    # not equal!
    test.gls <- residuals2(e.gls, adjust.residuals = FALSE)
    test.lme <- residuals2(e.lme, adjust.residuals = FALSE)
    expect_equal(unname(test.gls),unname(test.lme))
    
    test.lvm <- residuals2(e.lvm, adjust.residuals = FALSE)
    expect_equal(test.lvm,test.gls)

    GS.lme <- as.double(residuals(e.lme, type = "response", level = 0))
    GS.gls <- as.double(residuals(e.gls))
    expect_equal(GS.lme,GS.gls)
    
    expect_equal(GS.lme,as.double(test.gls))
    expect_equal(GS.lme,as.double(residuals(e.lvm)))
})

## ** versus lvm
m <- lvm(c(Y1~1*eta1,Y2~1*eta1,Y3~1*eta1,eta1~beta*G1,
           Z1~1*eta2,Z2~1*eta2,Z3~1*eta2,eta2~beta*G2)
         )
latent(m) <- ~eta1+eta2
d <- sim(m, 5e1)
e.lvm <- estimate(m,d)

e.lvm2 <- e.lvm
e.lvm2$prepareScore2 <- prepareScore2(lava::Model(e.lvm2), data = d)

test_that("equivalence residuals2.lvm residuals.lvm", {
    test <- residuals2(e.lvm, adjust.residuals = FALSE)
    test2 <- residuals2(e.lvm2, adjust.residuals = FALSE)    
    GS <- residuals(e.lvm)
    expect_equal(GS,test)
    expect_equal(GS,test2)
})



##----------------------------------------------------------------------
### test-residuals.R ends here

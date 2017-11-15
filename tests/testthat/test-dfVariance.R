### test-dfVariance.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 20 2017 (10:22) 
## Version: 
## last-updated: nov 15 2017 (17:58) 
##           By: Brice Ozenne
##     Update #: 55
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(testthat)
library(clubSandwich)
library(nlme)
library(lme4)
library(lmerTest)

context("dfVariance")
n <- 5e1

## * linear regression
mSim <- lvm(Y1~X1+X2+X3,Y2~X2,Y3~1)
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
d <- as.data.table(sim(mSim,n))

## ** t test
## formula:
## df = \frac{ 2 * s_pool^2 }{ var(s_pool^2) }
##    = \frac{ ( s_X^2/m + s_Y^2/n )^2}{( s_X^4/(m(m-1)) + s_Y^4/(n(n-1)))}

## using the t test function
e.ttest <- t.test(d$X1,d$X2)
e.ttest$parameter

## by hand
sX1 <- var(d$X1)/n
sX2 <- var(d$X2)/n
df <- (sX1+sX2)^2/(sX1^2/(n-1) + sX2^2/(n-1))

df-e.ttest$parameter

## ** lm
e.lvm <- estimate(lvm(Y1~X1+X2),d)
e.lm <- lm(Y1~X1+X2,d)

### *** clubSandwich
cS.vcov <- vcovCR(e.lm, type = "CR2", cluster = d$Id)
cS.df <- coef_test(e.lm, vcov = cS.vcov, test = "Satterthwaite", cluster = 1:NROW(d))
## cS.df$df is very suspect: should be the same for all coefficient and close to n-p

### *** dfVariance
test_that("linear regression: df",{

    df.adj <- dfVariance(e.lvm,
                         robust = FALSE, adjust.residuals = FALSE)

    keep.coef <- c("Y1","Y1~X1","Y1~X2")
    GS <- rep(NROW(d),length(keep.coef))
    expect_equal(as.double(df.adj[keep.coef]),GS)

    df.adj <- dfVariance(e.lvm,
                         robust = TRUE, adjust.residuals = FALSE)
    df.adj

})

## * mixed model
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G+Gender))
latent(mSim) <- ~eta
categorical(mSim, labels = c("M","F")) <- ~Gender
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- as.data.table(sim(mSim,n,latent = FALSE))
setkey(dW, "Id")
dL <- melt(dW,id.vars = c("G","Id","Gender"), variable.name = "time")
setkey(dL, "Id")

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

## ** clubSandwich - bug
expect_equal(logLik(e.lmer),logLik(e.lme))
coef_test(e.lme, vcov = "CR0", test = "Satterthwaite", cluster = dL$Id)
## strange that same type of coef have very different degrees of freedom

## ** lmerTest - ok
summary(e.lmer, ddf = "Satterthwaite")$coef

## ** lava
expect_equal(as.double(logLik(e.lmer)),as.double(logLik(e.lvm)))

test_that("mixed mode: df",{
    iidX <- iid2(e.lvm)
    iidY <- iid2(e.lme)
    range(iidX[,c(1:5,7,6)]-iidY)
    
    V1 <- solve(information(e.lvm, p = pars(e.lvm)+1))
    V2 <- attr(residuals2(e.lme, p = .coef2(e.lme)+1,
                          adjust.residuals = FALSE, return.vcov.param = TRUE),
               "vcov.param")
    V1-V2
    
    df.adj <- dfVariance(e.lvm,
                         robust = FALSE, adjust.residuals = FALSE)

    expect_equal(as.double(summary(e.lmer, ddf = "Satterthwaite")$coef[,"df"]),
                 as.double(df.adj[1:5]))

    residuals2(e.gls)
    df.lme <- dfVariance(e.lme,
                         robust = FALSE, adjust.residuals = FALSE)
    
    
    vcov.param <- attr(residuals2(e.lme, p = .coef2(e.lme)+1,
                                  adjust.residuals = FALSE, return.vcov.param = TRUE),
                       "vcov.param")
    vcov.param
    vcov(e.lme)
    

})

#----------------------------------------------------------------------
### test-dfVariance.R ends here


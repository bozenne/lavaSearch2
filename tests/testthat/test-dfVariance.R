### test-dfVariance.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 20 2017 (10:22) 
## Version: 
## last-updated: nov 16 2017 (16:52) 
##           By: Brice Ozenne
##     Update #: 70
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
## ** t test
## formula:
## df = \frac{ 2 * s_pool^2 }{ var(s_pool^2) }
##    = \frac{ ( s_X^2/m + s_Y^2/n )^2}{( s_X^4/(m(m-1)) + s_Y^4/(n(n-1)))}

## using the t test function
e.ttest <- t.test(dW$Y1,dW$Y2)
e.ttest$parameter

## by hand
sX1 <- var(dW$Y1)/n
sX2 <- var(dW$Y2)/n
df <- (sX1+sX2)^2/(sX1^2/(n-1) + sX2^2/(n-1))

df-e.ttest$parameter

## ** lm
e.lvm <- estimate(lvm(Y1~X1+X2), data = dW)
e.lm <- lm(Y1~X1+X2, data = dW)

### *** clubSandwich
cS.vcov <- vcovCR(e.lm, type = "CR2", cluster = d$Id)
cS.df <- coef_test(e.lm, vcov = cS.vcov, test = "Satterthwaite", cluster = 1:NROW(d))
cS.df
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

## *** clubSandwich - bug
expect_equal(logLik(e.lmer),logLik(e.lme))
coef_test(e.lme, vcov = "CR0", test = "Satterthwaite", cluster = dL$Id)
## strange that same type of coef have very different degrees of freedom

## *** lmerTest - ok
summary(e.lmer, ddf = "Satterthwaite")$coef
 
## *** lava - ok
expect_equal(as.double(logLik(e.lmer)),as.double(logLik(e.lvm)))

test_that("mixed mode: df",{
    GS <- summary(e.lmer, ddf = "Satterthwaite")$coef[,"df"]
    
    df.adj.lvm <- dfVariance(e.lvm,
                             robust = FALSE, adjust.residuals = FALSE)
    expect_equal(as.double(GS),
                 as.double(df.adj.lvm[1:5]))

    df.adj.lme <- dfVariance(e.lme,
                             robust = FALSE, adjust.residuals = FALSE)
    expect_equal(GS, df.adj.lme[names(GS)])

    df.adj.gls <- dfVariance(e.gls,
                             robust = FALSE, adjust.residuals = FALSE)
    expect_equal(GS, df.adj.gls[names(GS)], tol = 1e-4)
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

logLik(e.lvm)
logLik(e.lme)
logLik(e.gls)

test_that("mixed mode: df",{
    df.adj.lvm <- dfVariance(e.lvm, fast = FALSE,
                             robust = FALSE, adjust.residuals = FALSE)
    df.adj.lvm

    ## df.adj.lme <- dfVariance(e.lme,
    ##                          robust = FALSE, adjust.residuals = FALSE)

    df.adj.gls <- dfVariance(e.gls,
                             robust = FALSE, adjust.residuals = FALSE)
    df.adj.gls
})


#----------------------------------------------------------------------
### test-dfVariance.R ends here


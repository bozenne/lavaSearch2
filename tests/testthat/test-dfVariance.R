### test-dfVariance.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 20 2017 (10:22) 
## Version: 
## last-updated: dec 13 2017 (16:19) 
##           By: Brice Ozenne
##     Update #: 85
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


df.lvm <- dfVariance(e.lvm,
                     robust = FALSE, adjust.residuals = TRUE)

e.lm <- lm(Y1~X1+X2, data = dW)


### *** clubSandwich
cS.vcov <- vcovCR(e.lm, type = "CR0", cluster = dW$Id)
cS.df <- coef_test(e.lm, vcov = cS.vcov, test = "Satterthwaite", cluster = 1:NROW(dW))
cS.df
## cS.df$df is very suspect: should be the same for all coefficient and close to n-p

### *** dfVariance
test_that("linear regression: df",{

    df.lvm <- dfVariance(e.lvm, adjust.residuals = FALSE)

    
    n.param <- length(df.lvm)
    GS <- c(rep(NROW(dW),n.param-1), NROW(dW)/4)
    expect_equal(as.double(df.lvm),GS)
    
    df.adj.lvm <- dfVariance(e.lvm, robust = FALSE, adjust.residuals = TRUE)
    df.adj.lvm
    coef(e.lvm)
        
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

library(pbkrtest)

## *** lava - ok
expect_equal(as.double(logLik(e.lmer)),as.double(logLik(e.lvm)))

test_that("mixed mode: df",{
    GS <- summary(e.lmer, ddf = "Satterthwaite")$coef[,"df"]
    
    df.lvm <- dfVariance(e.lvm, adjust.residuals = FALSE)
    expect_equal(as.double(GS),
                 as.double(df.lvm[1:5]))

    df.lme <- dfVariance(e.lme, adjust.residuals = FALSE)
    expect_equal(GS, df.lme[names(GS)])

    df.gls <- dfVariance(e.gls, adjust.residuals = FALSE)
    expect_equal(GS, df.gls[names(GS)], tol = 1e-4)


    df.adj.lvm <- dfVariance(e.lvm, robust = FALSE, adjust.residuals = TRUE, fix.mean = TRUE)
    df.adj.lvm

    df.adj.lme <- dfVariance(e.lme, robust = FALSE, adjust.residuals = TRUE)
    df.adj.lme

    df.adj.gls <- dfVariance(e.gls, robust = FALSE, adjust.residuals = TRUE)
    df.adj.gls
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


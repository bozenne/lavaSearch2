### test-df.residuals.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 20 2017 (10:22) 
## Version: 
## last-updated: okt 24 2017 (16:50) 
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
library(clubSandwich)
library(nlme)
library(lme4)
library(lmerTest)

context("df.residuals")
n <- 5e1

## * linear regression
mSim <- lvm(Y1~X1+X2+X3,Y2~X2,Y3~1)
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
d <- as.data.table(sim(mSim,n))

## ** t test
## formula:
## df = \frac{ 2 * s_pool^2 }{ var(s_pool^2) }
##    = \frac{ ( s_X^2/m + s_Y^2/n )^2}{( s_X^4/(m(m-1)) + s_Y^4/(n(n-1)))}
e.ttest <- t.test(d$X1,d$X2)
e.ttest$parameter


sX1 <- var(d$X1)/n
sX2 <- var(d$X2)/n
df <- (sX1+sX2)^2/(sX1^2/(n-1) + sX2^2/(n-1))


dtL <- melt(d[,.(X1,X2,Id)], id.vars = "Id")
e.gls <- gls(value ~ 0+variable, weight = varIdent(form =~ 1|variable), data = dtL, method = "ML")
cS.df <- coef_test(e.gls, vcov = "CR2", test = "Satterthwaite", cluster = dtL$Id) # ok
cS.df <- coef_test(e.gls, vcov = "CR2", test = "Satterthwaite", cluster = 1:NROW(dtL)) # not ok

  epsilon = as.double(iIH %*% epsilon.tempo) 

## ** lm
e.lvm <- estimate(lvm(Y1~X1+X2),d)
e.lm <- lm(Y1~X1+X2,d)

cS.vcov <- vcovCR(e.lm, type = "CR2", cluster = d$Id)
cS.df <- coef_test(e.lm, vcov = cS.vcov, test = "Satterthwaite", cluster = 1:NROW(d))
## cS.df$df is very suspect: should be the same for all coefficient and close to n-p


test_that("linear regression: df",{

    s.lvm <- score2(e.lvm, adjust.residuals = TRUE)
    expect_equal(df.residual(e.lvm, adjust.residuals = TRUE),df.residual(e.lm))

    
    expect_equal(attr(s.lvm,"df"),df.residual(e.lm))

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

m <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G+Gender))
e.lvm <- estimate(m, dW)

e.lmer <- lmer(value ~ time + G + Gender + (1|Id), data = dL, REML = FALSE)
summary(e.lmer, ddf = "Satterthwaite")

e.lme <- lme(value ~ time + G + Gender, random = ~ 1|Id, data = dL, method = "ML")
expect_equal(logLik(e.lmer),logLik(e.lme))
coef_test(e.lme, vcov = "CR2", test = "Satterthwaite", cluster = dL$Id)
## strange that same type of coef have very different degrees of freedom

## * linear regressions
m.sim <- lvm(Y1~X1+X2+X3,Y2~X2,Y3~1)
d <- sim(m.sim,n)

e.lvm <- estimate(m.sim,d)

test_that("linear regressions: df",{

    s.lvm <- score2(e.lvm, adjust.residuals = FALSE)
    expect_equal(df.residual(e.lvm),attr(s.lvm,"df"))
    expect_equal(df.residual(e.lvm) ,)
    
})



e1 <- estimate(lvm(Y1~X1,Y2~X2+X3,Y3~1),d)
s1 <- score2(e1, adjust.residuals = TRUE)
expect_equal(s1[,c("Y1","Y1~X1","Y1~~Y1")],
                 s0[,c("Y1","Y1~X1","Y1~~Y1")])
    expect_equal(df.residual(e0),attr(s0,"df"))
    ## expect_equal(df.residual(e1),attr(s1,"df"))

    ## df.residual(e1)
    
##    expect_equal(attr(s1, "df"), NROW(d))


#----------------------------------------------------------------------
### test-df.residuals.R ends here

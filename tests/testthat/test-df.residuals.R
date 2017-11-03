### test-df.residuals.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 20 2017 (10:22) 
## Version: 
## last-updated: okt 27 2017 (14:37) 
##           By: Brice Ozenne
##     Update #: 29
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

lava.options(Dmethod = "Richardson")
iid.lava <- iid(e.lvm)
iid.manual <- score(e.lvm,indiv = TRUE) %*% vcov(e.lvm)
iid.lava-iid.manual
head(score(e.lvm,indiv = TRUE))
head(iid.lava)
crossprod(iid.lava)

cS.vcov <- vcovCR(e.lm, type = "CR2", cluster = d$Id)
cS.df <- coef_test(e.lm, vcov = cS.vcov, test = "Satterthwaite", cluster = 1:NROW(d))
## cS.df$df is very suspect: should be the same for all coefficient and close to n-p

score2(e.lvm)

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

m <- lvm(c(Y1[mu1:sigma]~1*eta,
           Y2[mu2:sigma]~1*eta,
           Y3[mu3:sigma]~1*eta,
           eta~G+Gender))
e.lvm <- estimate(m, dW)

## ** clubSandwich - bug
e.lme <- lme(value ~ time + G + Gender, random = ~ 1|Id, data = dL, method = "ML")
expect_equal(logLik(e.lmer),logLik(e.lme))
coef_test(e.lme, vcov = "CR0", test = "Satterthwaite", cluster = dL$Id)
## strange that same type of coef have very different degrees of freedom

## ** lmerTest - ok
e.lmer <- lmer(value ~ time + G + Gender + (1|Id),
               data = dL, REML = FALSE)
summary(e.lmer, ddf = "Satterthwaite")$coef


e.lme <- lme(value ~ time + G + Gender,
             random =~ 1|Id,
             weight = varIdent(form = ~1|Gender),
             data = dL, method = "ML")
class(e.lme$modelStruct$varStruct)
iid2(e.lme)
sqrt(diag(crossprod(iid2(e.lme, adjust.residuals = FALSE, return.df = FALSE))))
vcov(e.lme)
##             Estimate Std. Error       df t value Pr(>|t|)    
## (Intercept)  -0.5164     0.2496  80.3100  -2.068   0.0418 *  
## timeY2        0.3653     0.2063 100.0000   1.771   0.0796 .  
## timeY3        0.1273     0.2063 100.0000   0.617   0.5384    
## G             0.8964     0.1562  50.0000   5.738 5.58e-07 ***
## GenderF       1.4214     0.2891  50.0000   4.917 9.91e-06 ***

diag(rep(1,length(rho$fixEffs)))

ls.rho <- lmerTest:::rhoInit(list(), e.lmer, FALSE)
ls.rho$A <- lmerTest:::calcApvar(ls.rho)
n.param <- length(ls.rho$fixEffs)
L <- diag(rep(1,n.param))

iid2(e.lvm, adjust.residuals = FALSE)

score(e.lvm,indiv = TRUE) - score2(e.lvm, adjust.residuals = FALSE, Dmethod = "simple")

names(e.lvm)
names(e.lvm$model)
e.lvm$model$covpar
e.lvm$model$index
e.lvm$model$parpos

coefType(e.lvm)
coefType(m)

coef(e.lvm,level=9)

coef(e.lvm)
vss <- lmerTest:::vcovLThetaL(e.lmer)
vss(t(L[1,]), c(ls.rho$thopt,ls.rho$sigma))

undebug(vss)

logLik(e.lvm)
logLik(e.lmer)
logLik(e.lme)

lapply(1:n.param, function(iP){ # iP <- 1

    iL <- L[iP,]
    calcSatterth1DF2(ls.rho, L = iL, isF = FALSE)

    debug(vss)
    vss <- lmerTest:::vcovLThetaL(ls.rho$model)

    fct.obj <- function(x){
        print(x)
        vss(t(iL), x)
    }
    fct.obj <- function(x){ # x <- p.obj
        iid.tempo <- iid2(e.lme, p = x, adjust.residuals = FALSE, return.df = FALSE)
        vec.sd <- sqrt(diag(crossprod()))
        vec.sd %*% iL
    }
    p.obj <- c(ls.rho$thopt, ls.rho$sigma)    
    g <- as.double(numDeriv::jacobian(func = fct.obj, x = p.obj, method = "Richardson"))
        
    denom <- t(g) %*% ls.rho$A %*% g
    varcor <- fct.obj(p.obj)
    
    2 * (varcor)^2/denom

    
    result[, 2] <- (L %*% ls.rho$fixEffs)/sqrt(varcor)
    
    result[, 3] <- 2 * (1 - pt(abs(result[, 2]), df = result[, 
        1]))
    result[, 4] <- sqrt(varcor)
    
})

ls.rho$sigma

sigma(e.lmer)
e.lmer@theta


str(e.lmer)

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

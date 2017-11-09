### test-df.residuals.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 20 2017 (10:22) 
## Version: 
## last-updated: nov  9 2017 (18:19) 
##           By: Brice Ozenne
##     Update #: 35
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
logLik(e.lmer)

## *** try to replicate lmerTest results
e.lme <- lme(value ~ time + G + Gender,
             random =~ 1|Id,
             data = dL, method = "ML")
logLik(e.lme)

m.lvm <- lvm(Y1[mu1:sigma]~1*eta,Y2[mu2:sigma]~1*eta,Y3[mu3:sigma]~1*eta,eta~G+Gender)
latent(m.lvm) <- ~eta
e.lvm <- estimate(m.lvm, data = dW)
logLik(e.lvm)

solve(information(e.lvm))-attr(residuals2(e.lme, return.vcov.param = TRUE),"vcov.param")
solve(information(e.lvm))-attr(residuals2(e.lvm, return.vcov.param = TRUE),"vcov.param")



##              (Intercept)        timeY2        timeY3             G       GenderF     corCoef1       sigma2
## (Intercept)  0.062322690 -2.127369e-02 -2.127369e-02  2.764622e-03 -4.842052e-02  0.000000000  0.000000000
## timeY2      -0.021273689  4.254738e-02  2.127369e-02  6.554524e-19 -6.372999e-18  0.000000000  0.000000000
## timeY3      -0.021273689  2.127369e-02  4.254738e-02  4.574516e-19 -3.214343e-18  0.000000000  0.000000000
## G            0.002764622  6.554524e-19  4.574516e-19  2.441123e-02 -5.239548e-03  0.000000000  0.000000000
## GenderF     -0.048420521 -6.372999e-18 -3.214343e-18 -5.239548e-03  8.358517e-02  0.000000000  0.000000000
## corCoef1     0.000000000  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.042864621 -0.007542831
## sigma2       0.000000000  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -0.007542831  0.022628493

param <- coef(e.lvm)
n.param <- length(param)
L <- vector("numeric", length = n.param)
L[4] <- 1
fct.obj <- function(x){ # x <- p.obj
    newSigma <- solve(information(e.lvm, p = x))
    return(rbind(L) %*% newSigma %*% cbind(L))
    ## iid.tempo <- iid2(e.lme, p = x, adjust.residuals = FALSE, return.df = FALSE)
    ## vec.sd <- sqrt(diag(crossprod()))
    ## vec.sd %*% iL
}
fct.obj(x = param)

g <- as.double(numDeriv::jacobian(func = fct.obj, x = param, method = "Richardson"))     
denom <- t(g) %*% vcov(e.lvm)  %*% g
2*fct.obj(x = param)^2 / (t(g) %*% vcov(e.lvm)  %*% g)
calcDF_lmerTest(e.lmer)



    

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


## * Robust vcov

## ** linear regression
mSim <- lvm(Y~X1+X2+X3)
transform(mSim,Id~Y) <- function(x){1:NROW(x)}
set.seed(10)
d <- sim(mSim,n)

m <- lvm(Y~X1+X2+X3)
e <- estimate(m,d)
eR <- estimate(m,d,robust = TRUE, cluster = "Id")

test_that("linear regression (at ML)",{

    test <- solve(vcov(e))
    X <- model.matrix(lm(formula(m)[[1]], data=d))    
    GS <- Matrix::bdiag(crossprod(X)/coef(e)["Y~~Y"],n/2*(coef(e)["Y~~Y"])^(-2))
    expect_equal(as.double(test), as.double(GS))

    test <- crossprod(iid(e, return.df = FALSE))
    
    GS <- rmAttr(vcov(e),c("det","pseudo","minSV"))
    expect_equal(unname(test), unname(GS))    
})


p <- length(coef(e))
fSigma <- function(p){
    vcovB <- p["Y~~Y"]*solve(crossprod(X))
    vcovS <- 2/n*(p["Y~~Y"])^2   
    
    return(as.double(Matrix::bdiag(vcovB,vcovS)))
}
matrix(fSigma(coef(e)),p,p) - vcov(e)

iGrad <- numDeriv::jacobian(func = fSigma, x = coef(e))
iGrad[,5]
solve(crossprod(X))
(4/n*(coef(e)["Y~~Y"]))

(2/n*(coef(e)["Y~~Y"])^2) * solve(crossprod(X))[1] * solve(crossprod(X))[1]

2 / (solve(crossprod(X))[1] * solve(crossprod(X))[1])

iGrad %*% vcov(e) %*% iGrad

2/n*(p["Y~~Y"])^2   

#----------------------------------------------------------------------
### test-df.residuals.R ends here

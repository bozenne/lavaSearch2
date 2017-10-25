### test-score2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 13 2017 (11:28) 
## Version: 
## last-updated: okt 25 2017 (12:06) 
##           By: Brice Ozenne
##     Update #: 117
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
n <- 3e1

## * score for nlme models
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
latent(mSim) <- ~eta
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- as.data.table(sim(mSim,n,latent = FALSE))
setkey(dW, "Id")
dL <- melt(dW,id.vars = c("G","Id"), variable.name = "time")
setkey(dL, "Id")
dL$Z1 <- rnorm(NROW(dL))

## ** gls models (Heteroscedasticity)
m <- lvm(c(Y1~G,Y2~G,Y3~G))
e.lvm <- estimate(m, dW)

e.gls <- gls(value ~ 0+time + time:G,
             weight = varIdent(form = ~ 1|time),
             data = dL, method = "ML")

keep.cols <- c("Y1","Y2","Y3","Y1~G","Y2~G","Y3~G")

test_that("lme equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.gls)))
})

test_that("score2 equivalent to score", {
    score.gls <- score2(e.gls, cluster = "Id", adjust.residuals = FALSE, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = FALSE, indiv = TRUE)

    expect_equal(unname(score.gls),unname(score.lvm[,keep.cols]), tol = 1e-5)
    
    score.gls <- score2(e.gls, cluster = "Id", adjust.residuals = TRUE, power = 1, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = TRUE, power = 1, indiv = TRUE)

    score.lvm[,keep.cols]/score.gls
    vcov(e.gls)/vcov(e.lvm)[keep.cols,keep.cols]

    score.gls <- score2(e.gls, cluster = "Id", adjust.residuals = TRUE, power = 0.5, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = TRUE, power = 0.5, indiv = TRUE)

    score.lvm[,keep.cols]/score.gls
})

## ** lme
m <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
e.lvm <- estimate(m, dW)

e.lme <- lme(value ~ time + G,
             random =~1| Id,
             weight = varIdent(form = ~ 1|time),
             data = dL, method = "ML")

test_that("lme equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
  
    coef.lme <- c(fixef(e.lme), sigma(e.lme)^2, as.numeric(getVarCov(e.lme)),
    (sigma(e.lme)*coef(e.lme$modelStruct$varStruct, uncons = FALSE, allCoef = FALSE))^2)
    coef.lvm <- coef(e.lvm)

    expect_equal(as.double(coef.lme),as.double(coef.lvm), tol = 1e-5)
})

test_that("score2 equivalent to score", {
    score.lme <- score2(e.lme, adjust.residuals = FALSE, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = FALSE, indiv = TRUE)

    expect_equal(unname(score.lme),unname(score.lvm[,c("eta","Y2","Y3","eta~G")]), tol = 1e-5)

    score.lme <- score2(e.lme, adjust.residuals = TRUE, indiv = TRUE, power = 1)
    score.lvm <- score2(e.lvm, adjust.residuals = TRUE, indiv = TRUE, power = 1)

    expect_equal(unname(score.lme),unname(score.lvm[,c("eta","Y2","Y3","eta~G")]), tol = 1e-5)
    
    score.lme <- score2(e.lme, adjust.residuals = TRUE, indiv = TRUE, power = 0.5)
    score.lvm <- score2(e.lvm, adjust.residuals = TRUE, indiv = TRUE, power = 0.5)

    expect_equal(unname(score.lme),unname(score.lvm[,c("eta","Y2","Y3","eta~G")]), tol = 1e-5)
})

## ** gls vs. lme models 
e.gls <- gls(value ~ time*G,
             correlation = corCompSymm(form =~ 1 | Id),
             data = dL, method = "ML")

e.lme <- lme(value ~ time*G,
             random =~1| Id,
             data = dL, method = "ML")

test_that("lme equivalent to gls", {
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lme)))
  
    tau.lme <- as.numeric(getVarCov(e.lme))
    sigma2.lme <- as.numeric(getVarCov(e.lme)) + sigma(e.lme)^2
    coef.lme <- c(fixef(e.lme), sigma2 = sigma2.lme, tau = tau.lme)

    tau.gls <- getVarCov(e.gls)[1,2]
    sigma2.gls <- getVarCov(e.gls)[1,1]
    coef.gls <- c(coef(e.gls), sigma2 =  sigma2.gls, tau = tau.gls)
    
    expect_equal(as.double(coef.lme),as.double(coef.gls), tol = 1e-7)

    score.lme <- score2(e.lme, adjust.residuals = FALSE, indiv = TRUE)
    score.gls <- score2(e.gls, adjust.residuals = FALSE, indiv = TRUE)

    expect_equal(unname(score.lme),unname(score.gls), tol = 1e-5)

    score.lme <- score2(e.lme, adjust.residuals = TRUE, indiv = TRUE, power = 1)
    score.gls <- score2(e.gls, adjust.residuals = TRUE, indiv = TRUE, power = 1)
        
    expect_equal(unname(score.lme),unname(score.gls), tol = 1e-5)
    
    score.lme <- score2(e.lme, adjust.residuals = TRUE, indiv = TRUE, power = 0.5)
    score.gls <- score2(e.gls, adjust.residuals = TRUE, indiv = TRUE, power = 0.5)

    expect_equal(unname(score.lme),unname(score.gls), tol = 1e-5)
})


## * not-adjusted score
## ** linear regression
m <- lvm(Y~X1+X2+X3)
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)

test_that("linear regression",{
    ## at ML
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("linear regression: constrains",{
    m <- lvm(Y~X1+1*X2)
    e <- estimate(m, d)

    expect_equal(score2(e, adjust.residuals = FALSE),
                 score(e, indiv = TRUE))
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
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

## *** with covariance links
m <- lvm(c(Y1~X1,Y2~X2,Y3~X3+X1))
covariance(m) <- Y1~Y2
set.seed(10)
d <- sim(m,n)

e <- estimate(m,d)
param <- coef(e)

test_that("multiple linear regression (covariance link)",{
        ## at ML
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("linear regressions: constrains",{
    m <- lvm(Y1~X1+1*X2,Y2~2*X3+2*X1,Y3~X2)
    e <- estimate(m, d)

    expect_equal(score2(e, adjust.residuals = FALSE),
                 score(e, indiv = TRUE))
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
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))
})

test_that("factor model: constrains",{
    m <- lvm(Y1~1*eta+1*X2,Y2~1*eta,Y3~1*eta)
    e <- estimate(m, d)

    expect_equal(score2(e, adjust.residuals = FALSE),
                 score(e, indiv = TRUE))
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
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
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
    test <- score2(e, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, indiv=TRUE)
    expect_equal(test, GS)
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
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
    
    test <- score2(e, indiv = FALSE, adjust.residuals = FALSE)
    GS <- score(e, indiv=FALSE)
    expect_equal(as.double(test),as.double(GS))

    ## not at ML
    test <- score2(e, p = param+1, indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+1, indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+0.1*(1:length(param)), indiv=TRUE, adjust.residuals = FALSE)
    GS <- score(e, p = param+0.1*(1:length(param)), indiv=TRUE)
    expect_equal(test, GS)

    test <- score2(e, p = param+1, indiv = FALSE, adjust.residuals = FALSE)
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

## ** multigroup model
m.sim <- lvm(Y~X1+X2,G~1)
categorical(m.sim,K=2,label=c("a","b")) <- ~G+X2
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)
e <- estimate(list(m.sim,m.sim),data = split(d,d$G))

expect_error(score2(e))

## ** model with binary endogenous variables
library(lava.tobit)
m.sim <- lvm(Y~X1)
categorical(m.sim,K=2,labels = c("a","b")) <- ~Y
set.seed(10)
d <- sim(m.sim,n,latent=FALSE)
e <- estimate(lvm(Y~X1),data = d)

expect_error(score2(e))

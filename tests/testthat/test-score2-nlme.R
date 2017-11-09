### test-score2-nlme.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: nov  6 2017 (11:40) 
## Version: 
## last-updated: nov  9 2017 (16:27) 
##           By: Brice Ozenne
##     Update #: 33
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


library(testthat)
library(nlme)

context("iid2-nlme")

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

## ** gls models (Homoschedasticity)
m <- lvm(c(Y1~G))
e.lvm <- estimate(m, dW)

e.gls <- gls(Y1 ~ G, data = dW, method = "ML")
allCoef.gls <- c(coef(e.gls),sigma2 = sigma(e.gls)^2)

test_that("gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.gls)))
})

test_that("score2.gls equivalent to score.lvm", {
    score.gls <- score2(e.gls, cluster = "Id", adjust.residuals = FALSE, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = FALSE, indiv = TRUE)

    expect_equal(unname(score.gls),unname(score.lvm))    
})

## ** gls models (Heteroscedasticity)
m <- lvm(c(Y1~G,Y2~G,Y3~G))
e.lvm <- estimate(m, dW)

e.gls <- gls(value ~ 0+time + time:G,
             weight = varIdent(form = ~ 1|time),
             data = dL, method = "ML")
allCoef.gls <- c(coef(e.gls),sigma2 = sigma(e.gls)^2,coef(e.gls$modelStruct$varStruct, unconstrained = FALSE, allCoef = FALSE))

test_that("gls equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.gls)))
})

test_that("score2.gls equivalent to score.lvm", {

    score.gls <- score2(e.gls, cluster = "Id", adjust.residuals = FALSE, indiv = TRUE)
    score.gls.p <- score2(e.gls, p = allCoef.gls, cluster = "Id", adjust.residuals = FALSE, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = FALSE, indiv = TRUE)
    
    expect_equal(unname(score.gls[,1:6]),unname(score.lvm[,1:6]))
    expect_true(all(abs(colSums(score.gls))<1e-7))
    expect_equal(score.gls,score.gls.p)
})

test_that("score2.gls equivalent to score.lvm (adjust residuals)", {
    score.gls <- score2(e.gls, cluster = "Id", adjust.residuals = TRUE, power = 1, indiv = TRUE)
    score.gls.p <- score2(e.gls, p = allCoef.gls, cluster = "Id", adjust.residuals = TRUE, power = 1, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = TRUE, power = 1, indiv = TRUE)

    expect_equal(unname(score.gls[,1:6]),unname(score.lvm[,1:6]))
    expect_equal(score.gls,score.gls.p)
    
    score.gls <- score2(e.gls, cluster = "Id", adjust.residuals = TRUE, power = 0.5, indiv = TRUE)
    score.gls.p <- score2(e.gls, p = allCoef.gls, cluster = "Id", adjust.residuals = TRUE, power = 0.5, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = TRUE, power = 0.5, indiv = TRUE)

    expect_equal(unname(score.gls[,1:6]),unname(score.lvm[,1:6]))
    expect_equal(score.gls,score.gls.p)
})

## ** lme
m <- lvm(c(Y1[mu1:sigma2]~1*eta,Y2[mu2:sigma2]~1*eta,Y3[mu3:sigma2]~1*eta,eta~G))
e.lvm <- estimate(m, dW)

e.lme <- lme(value ~ time + G,
             random =~1| Id,
             data = dL, method = "ML")

allCoef.lme <- c(fixef(e.lme),sigma2 = sigma(e.lme)^2, corCoef1 = as.double(getVarCov(e.lme)))

test_that("lme equivalent to lvm", {
    expect_equal(as.double(logLik(e.lvm)), as.double(logLik(e.lme)))
  
    coef.lme <- c(fixef(e.lme), sigma(e.lme)^2, as.numeric(getVarCov(e.lme)),
    (sigma(e.lme)*coef(e.lme$modelStruct$varStruct, uncons = FALSE, allCoef = FALSE))^2)
    coef.lvm <- coef(e.lvm)

    expect_equal(as.double(coef.lme),as.double(coef.lvm))
})

test_that("score2.lme equivalent to score2.lvm", {
    score.lme <- score2(e.lme, adjust.residuals = FALSE, indiv = TRUE)
    score.lme.p <- score2(e.lme, p = allCoef.lme, adjust.residuals = FALSE, indiv = TRUE)
    score.lvm <- score2(e.lvm, adjust.residuals = FALSE, indiv = TRUE)
    # colSums(abs(score.lvm - score(e.lvm, indiv = TRUE)))

    expect_equal(unname(score.lme[,c(1:4,6,5)]),unname(score.lvm))
    expect_equal(score.lme,score.lme.p)
    
    score.lme <- score2(e.lme, adjust.residuals = TRUE, indiv = TRUE, power = 1)
    score.lme.p <- score2(e.lme, p = allCoef.lme, adjust.residuals = TRUE, indiv = TRUE, power = 1)
    score.lvm <- score2(e.lvm, adjust.residuals = TRUE, indiv = TRUE, power = 1)

    expect_equal(unname(score.lme[,c(1:4,6,5)]),unname(score.lvm))
    expect_equal(score.lme,score.lme.p)
    
    score.lme <- score2(e.lme, adjust.residuals = TRUE, indiv = TRUE, power = 0.5)
    score.lme.p <- score2(e.lme, p = allCoef.lme, adjust.residuals = TRUE, indiv = TRUE, power = 0.5)
    score.lvm <- score2(e.lvm, adjust.residuals = TRUE, indiv = TRUE, power = 0.5)

    expect_equal(unname(score.lme[,c(1:4,6,5)]),unname(score.lvm))
    expect_equal(score.lme,score.lme.p)
})

## ** gls vs. lme models 
e.gls <- gls(value ~ time*G,
             correlation = corCompSymm(form =~ 1 | Id),
#             weight = varIdent(form = ~ 1|time),
             data = dL, method = "ML")

e.lme <- lme(value ~ time*G,
             random =~1| Id,
#             weight = varIdent(form = ~ 1|time),
             data = dL, method = "ML")

name.coef <- names(fixef(e.lme))

test_that("lme equivalent to gls", {
    expect_equal(as.double(logLik(e.gls)), as.double(logLik(e.lme)))

    tau.lme <- as.numeric(getVarCov(e.lme))
    sigma2.lme <- as.numeric(getVarCov(e.lme)) + sigma(e.lme)^2
    coef.lme <- c(fixef(e.lme), sigma2 = sigma2.lme, tau = tau.lme)

    tau.gls <- getVarCov(e.gls)[1,2]
    sigma2.gls <- getVarCov(e.gls)[1,1]
    coef.gls <- c(coef(e.gls), sigma2 =  sigma2.gls, tau = tau.gls)
    
    expect_equal(as.double(coef.lme),as.double(coef.gls), tol = 1e-7)
})

test_that("score.lme equivalent to score.gls", {
    score.lme <- score2(e.lme, adjust.residuals = FALSE, indiv = TRUE)
    score.gls <- score2(e.gls, adjust.residuals = FALSE, indiv = TRUE)

    expect_equal(unname(score.lme[,name.coef]),unname(score.gls[,name.coef]), tol = 1e-7)
    
    score.lme <- score2(e.lme, adjust.residuals = TRUE, indiv = TRUE, power = 1)
    score.gls <- score2(e.gls, adjust.residuals = TRUE, indiv = TRUE, power = 1)

    expect_equal(unname(score.lme[,name.coef]),unname(score.gls[,name.coef]), tol = 1e-7)

    score.lme <- score2(e.lme, adjust.residuals = TRUE, indiv = TRUE, power = 0.5)
    score.gls <- score2(e.gls, adjust.residuals = TRUE, indiv = TRUE, power = 0.5)

    expect_equal(unname(score.lme[,name.coef]),unname(score.gls[,name.coef]), tol = 1e-7)
})



#----------------------------------------------------------------------
### test-score2-nlme.R ends here

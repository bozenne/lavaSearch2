### test1-sCorrect-ssc.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  7 2018 (12:08) 
## Version: 
## Last-Updated: jan 17 2020 (13:39) 
##           By: Brice Ozenne
##     Update #: 84
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
rm(list = ls())
if(FALSE){ ## already called in test-all.R
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))
library(nlme)
calcFactor <- function(object){
    return((object$dims$N - object$dims$p)/(object$dims$N - object$dims$p * (object$method == "REML")))
}

context("sCorrect (small sample correction)")

## * simulation
n <- 5e1
mSim <- lvm(c(Y1~eta1,Y2~eta1+X2,Y3~eta1+X1,
              Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(mSim) <- eta1~X1+Gender
latent(mSim) <- ~eta1+eta2
categorical(mSim, labels = c("Male","Female")) <- ~Gender
transform(mSim, Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
d <- lava::sim(mSim, n = n, latent = FALSE)
dL <- reshape2::melt(d, id.vars = c("Id","X1","X2","X3","Gender"),
                     measure.vars = c("Y1","Y2","Y3","Z1","Z2","Z3"))
dLred <- dL[dL$variable %in% c("Y1","Y2","Y3"),]
dLred$variable.factor <- as.factor(dLred$variable)

## * linear regression
## ** univariate
e.lm <- lm(Y1~X1+X2, data = d)
e.lvm <- estimate(lvm(Y1~X1+X2), data = d)

test_that("linear regression - residual correction equivalent to REML", {
    eSSC1.lm <- sCorrect(e.lm, ssc = "residuals")
    eSSC1.lvm <- sCorrect(e.lvm, ssc = "residuals")

    GS <- c(coef(e.lm), sigma(e.lm)^2)
    expect_equal(unname(eSSC1.lm$sCorrect$param),
                 unname(GS), tol = 1e-6)
    expect_equal(unname(eSSC1.lvm$sCorrect$param),
                 unname(GS), tol = 1e-6)

    GS <- vcov(e.lm)
    expect_equal(unname(eSSC1.lm$sCorrect$vcov.param[1:3,1:3]),
                 unname(GS), tol = 1e-6)
    expect_equal(unname(eSSC1.lvm$sCorrect$vcov.param[1:3,1:3]),
                 unname(GS), tol = 1e-6)
})

test_that("linear regression - Cox correction equivalent to REML", {
    eSSC2.lm <- sCorrect(e.lm, ssc = "Cox")
    eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")

    JJK <- array(0, dim = rep(p,3), dimnames = list(name.param,name.param,name.param))
    name.param <- c(names(coef(e.lm)),"sigma2")
    p <- length(name.param)
    X <- model.matrix(e.lm)

    JJK[name.param[1:3],name.param[1:3],"sigma2"] <- -crossprod(X)/sigma(e.lm)^4
    JJK[name.param[1:3],"sigma2",name.param[1:3]] <- -crossprod(X)/sigma(e.lm)^4
    JJK["sigma2",name.param[1:3],name.param[1:3]] <- crossprod(X)/sigma(e.lm)^4
    expect_equal(JJK, eSSC2.lm$sCorrect$ssc$JJK, tol = 1e-5)
    
    GS <- c(coef(e.lm), sigma(e.lm)^2)
    expect_equal(unname(eSSC2.lm$sCorrect$param),
                 unname(GS), tol = 1e-6)
    expect_equal(unname(eSSC2.lvm$sCorrect$param),
                 unname(GS), tol = 1e-6)

    GS <- vcov(e.lm)
    expect_equal(unname(eSSC2.lm$sCorrect$vcov.param[1:3,1:3]),
                 unname(GS), tol = 1e-6)
    expect_equal(unname(eSSC2.lvm$sCorrect$vcov.param[1:3,1:3]),
                 unname(GS), tol = 1e-6)
})


## ** multiple, no constrain
e.gls <- gls(value ~ -1 + variable + variable:X1,
             data = dLred,
             weight = varIdent(form = ~1|variable),
             method = "REML")
e.lvm <- estimate(lvm(Y1 ~ X1,
                      Y2 ~ X1,
                      Y3 ~ X1),
                  data = d)

test_that("multiple linear regression - residual correction equivalent to REML", {
    eSSC1.lvm <- sCorrect(e.lvm, ssc = "residuals")

    GS <- c(intervals(e.gls)$coef[,2],
            (c(Y1 = 1, intervals(e.gls)$varStruct[,2])*intervals(e.gls)$sigma[2])^2
            )
    
    expect_equal(unname(eSSC1.lvm$sCorrect$param[c("Y1","Y2","Y3","Y1~X1","Y2~X1","Y3~X1","Y1~~Y1","Y2~~Y2","Y3~~Y3")]),
                 unname(GS), tol = 1e-6)

    GS <- vcov(e.gls)
    expect_equal(unname(eSSC1.lvm$sCorrect$vcov.param[1:6,1:6]),
                 unname(GS), tol = 1e-6)
})

test_that("multiple linear regression - Cox correction equivalent to REML", {
    eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")

    GS <- c(intervals(e.gls)$coef[,2],
            (c(Y1 = 1, intervals(e.gls)$varStruct[,2])*intervals(e.gls)$sigma[2])^2
            )
    
    expect_equal(unname(eSSC2.lvm$sCorrect$param[c("Y1","Y2","Y3","Y1~X1","Y2~X1","Y3~X1","Y1~~Y1","Y2~~Y2","Y3~~Y3")]),
                 unname(GS), tol = 1e-6)

    GS <- vcov(e.gls)
    expect_equal(unname(eSSC2.lvm$sCorrect$vcov.param[1:6,1:6]),
                 unname(GS), tol = 1e-6)
})

## ** multiple, with constrains
e.gls0 <- gls(value ~ variable-1 + X1,
             data = dLred,
             weight = varIdent(form = ~1|variable),
             method = "ML")
e.gls <- gls(value ~ variable-1 + X1,
             data = dLred,
             weight = varIdent(form = ~1|variable),
             method = "REML")
e.lvm <- estimate(lvm(Y1[mu1:sigma1]~ beta1*X1,
                      Y2[mu2:sigma2]~ beta1*X1,
                      Y3[mu3:sigma3]~ beta1*X1),
                  data = d)
## logLik(e.gls)
## logLik(e.lvm)

test_that("multiple linear regression - residual correction equivalent to REML", {
    eSSC1.lvm <- sCorrect(e.lvm, ssc = "residuals")

    GS <- c(intervals(e.gls)$coef[,2],
            (c(Y1 = 1, intervals(e.gls)$varStruct[,2]) * intervals(e.gls)$sigma[2])^2)

    ## not precisely the same but better
    expect_equal(unname(eSSC1.lvm$sCorrect$param),
                 unname(GS), tol = 1e-3)
})

test_that("multiple linear regression - Cox correction equivalent to REML", {
    eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")

    GS <- c(intervals(e.gls)$coef[,2],
            (c(Y1 = 1, intervals(e.gls)$varStruct[,2]) * intervals(e.gls)$sigma[2])^2)

    ## not precisely the same but better
    expect_equal(unname(eSSC2.lvm$sCorrect$param),
                 unname(GS), tol = 1e-3)

})


## * mixed model
## ** CS
m <- lvm(c(Y1[0:sigma]~1*eta,
           Y2[0:sigma]~1*eta,
           Y3[0:sigma]~1*eta,
           eta~X1+X2))
latent(m) <- ~eta
e.lvm <- estimate(m, d)

e.lme <- nlme::lme(value ~ X1 + X2,
                   random =~1| Id,
                   data = dLred, method = "REML")

e.gls <- nlme::gls(value ~ X1 + X2,
                   correlation = corCompSymm(form = ~1| Id),
                   data = dLred, method = "REML")
 
test_that("mixed model (CS) - residual correction equivalent to REML", {
    eSSC1.lvm <- sCorrect(e.lvm, ssc = "residuals")

    GS <- c(intervals(e.lme)$fixed[,2],
            sigma2 = as.double(intervals(e.lme)$sigma[2])^2,
            tau = intervals(e.lme)$reStruct$Id[,2,]^2)

    ## not precisely the same but better
    expect_equal(unname(eSSC1.lvm$sCorrect$param),
                 unname(GS), tol = 1e-4)
    ## eSSC.lvm$sCorrect$param - GS
})

test_that("mixed model (CS) - Cox correction equivalent to REML", {
    eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")
    eSSC2.gls <- sCorrect(e.gls, ssc = "Cox")
    eSSC2.lme <- sCorrect(e.lme, ssc = "Cox")

    GS <- c(intervals(e.lme)$fixed[,2],
            sigma2 = as.double(intervals(e.lme)$sigma[2])^2,
            tau = intervals(e.lme)$reStruct$Id[,2,]^2)
    GS2 <- c(intervals(e.gls)$coef[,2],
            sigma2 = sigma(e.gls)^2,
            corCoef1 = intervals(e.gls)$corStruct[1,2])

    ## not precisely the same but better
    expect_equal(unname(eSSC2.lvm$sCorrect$param),
                 unname(GS), tol = 1e-4)
    expect_equal(unname(eSSC2.lme$sCorrect$param),
                 unname(GS), tol = 1e-4)
    expect_equal(unname(eSSC2.gls$sCorrect$param),
                 unname(GS2), tol = 1e-2)
    ## eSSC.lvm$sCorrect$param - GS
})

## ** CS with different variances
m <- lvm(c(Y1[0:sigma1]~1*eta,
           Y2[0:sigma2]~1*eta,
           Y3[0:sigma3]~1*eta,
           eta~X1+X2))
latent(m) <- ~eta
e.lvm <- estimate(m, d)

e.lme <- nlme::lme(value ~ X1 + X2,
                   random =~1| Id,
                   weights = varIdent(form =~ 1|variable),
                   data = dLred, method = "REML")

e.gls <- nlme::gls(value ~ X1 + X2,
                   correlation = corCompSymm(form = ~1| Id),
                   weights = varIdent(form =~ 1|variable),
                   data = dLred, method = "REML")

test_that("mixed model (CS,weight) - residual correction equivalent to REML", {
    eSSC1.lvm <- sCorrect(e.lvm, ssc = "residuals")

    GS <- c(intervals(e.lme)$fixed[,2],
            Y1 = as.double(intervals(e.lme)$sigma[2])^2,
            tau = intervals(e.lme)$reStruct$Id[,2,]^2,
            intervals(e.lme)$varStruct[,2]^2 * as.double(intervals(e.lme)$sigma[2])^2)

    ## not precisely the same but better
    expect_equal(unname(eSSC1.lvm$sCorrect$param),
                 unname(GS), tol = 5e-3)
    ## eSSC.lvm$sCorrect$param - GS
    ## coef(e.lvm) - GS
})

test_that("mixed model (CS,weight) - Cox correction equivalent to REML", {
    eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")
    eSSC2.gls <- sCorrect(e.gls, ssc = "Cox")
    eSSC2.lme <- sCorrect(e.lme, ssc = "Cox")

    GS <- c(intervals(e.lme)$fixed[,2],
            Y1 = as.double(intervals(e.lme)$sigma[2])^2,
            tau = intervals(e.lme)$reStruct$Id[,2,]^2,
            intervals(e.lme)$varStruct[,2]^2 * as.double(intervals(e.lme)$sigma[2])^2)

    GS2 <- c(intervals(e.lme)$fixed[,2],
             c(Y1 = 1,intervals(e.lme)$varStruct[,2])^2 * as.double(intervals(e.lme)$sigma[2])^2,
             tau = intervals(e.lme)$reStruct$Id[1,2])

    ## not precisely the same but better
    expect_equal(unname(eSSC2.lvm$sCorrect$param),
                 unname(GS), tol = 5e-3)
    eSSC2.lme$sCorrect$param - GS2
    ## eSSC.lvm$sCorrect$param - GS
})


## ** UN
m <- lvm(c(Y1~0+1*eta,
           Y2~0+1*eta,
           Y3~0+1*eta,
           eta~X1+X2))
covariance(m) <- Y1~Y2
covariance(m) <- Y1~Y3
e.lvm <- estimate(m, d)

e.gls <- nlme::gls(value ~ X1 + X2,
                   correlation = corSymm(form =~ 1| Id),
                   weight = varIdent(form =~ 1|variable),
                   data = dLred, method = "REML")

e.lme <- nlme::lme(value ~ X1 + X2,
                   random =~ 1|Id,
                   correlation = corSymm(),
                   weight = varIdent(form =~ 1|variable),
                   data = dLred, method = "REML")

test_that("mixed model (UN) - residual correction equivalent to REML", {
    eSSC1.lvm <- sCorrect(e.lvm, ssc = "residuals") 
    ## eSSC1.lvm$sCorrect$param - coef(e.lvm)
    
    gls_sigma2 <- as.double(intervals(e.gls)$sigma[2])^2
    gls_var <- c(Y1 = gls_sigma2, gls_sigma2 * intervals(e.gls)$varStruct[,2]^2)
    gls_tau <- as.double(sqrt(gls_var)["Y2"] * sqrt(gls_var)["Y3"] * intervals(e.gls)$corStruct[3,2])
    
    ## getVarCov2(eSSC1.lvm)
    ## eSSC1.lvm$sCorrect$param

    GS <- c(intervals(e.gls)$coef[,2],
            Y1 = as.double(gls_var["Y1"]) - gls_tau,
            "eta~~eta" = gls_tau,
            Y2 = as.double(gls_var["Y2"] - gls_tau),
            Y3 = as.double(gls_var["Y3"] - gls_tau),
            "Y1~~Y2" = as.double(sqrt(gls_var)["Y1"] * sqrt(gls_var)["Y2"] * intervals(e.gls)$corStruct[1,2] - gls_tau),
            "Y1~~Y3" = as.double(sqrt(gls_var)["Y1"] * sqrt(gls_var)["Y3"] * intervals(e.gls)$corStruct[2,2] - gls_tau)
            )

    ## not precisely the same but better
    expect_equal(unname(eSSC1.lvm$sCorrect$param),
                 unname(GS), tol = 5e-3)
    ## eSSC1.lvm$sCorrect$param - GS
    ## coef(e.lvm) - GS
})

test_that("mixed model (UN) - Cox correction equivalent to REML", {
    eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")

    ## eSSC2.lvm$sCorrect$param - coef(e.lvm)
    
    gls_sigma2 <- as.double(intervals(e.gls)$sigma[2])^2
    gls_var <- c(Y1 = gls_sigma2, gls_sigma2 * intervals(e.gls)$varStruct[,2]^2)
    gls_tau <- as.double(sqrt(gls_var)["Y2"] * sqrt(gls_var)["Y3"] * intervals(e.gls)$corStruct[3,2])

    GS <- c(intervals(e.gls)$coef[,2],
            Y1 = as.double(gls_var["Y1"]) - gls_tau,
            "eta~~eta" = gls_tau,
            Y2 = as.double(gls_var["Y2"] - gls_tau),
            Y3 = as.double(gls_var["Y3"] - gls_tau),
            "Y1~~Y2" = as.double(sqrt(gls_var)["Y1"] * sqrt(gls_var)["Y2"] * intervals(e.gls)$corStruct[1,2] - gls_tau),
            "Y1~~Y3" = as.double(sqrt(gls_var)["Y1"] * sqrt(gls_var)["Y3"] * intervals(e.gls)$corStruct[2,2] - gls_tau)
            )

    ## not precisely the same but better
    expect_equal(unname(eSSC2.lvm$sCorrect$param),
                 unname(GS), tol = 5e-3)
    ## eSSC2.lvm$sCorrect$param - GS
    ## coef(e.lvm) - GS
})

## * latent variable model
## ** factor model
m <- lvm(Y1~eta,
         Y2~eta+X2,
         Y3~eta,
         Z1~eta, Z1~~Y1,Z1~~Y2,
         eta~X1+X3)
e.lvm <- estimate(m, d)

## round(coef(estimate(m, sim(m,1e4))),1) ## truth

test_that("factor model - residuals correction", {
    eSSC1.lvm <- sCorrect(e.lvm, ssc = "residuals")
    ## coef(eSSC1.lvm) - coef(e.lvm)
    
    GS <- c("eta" = 0.23990945, "Y2" = 0.28804288, "Y3" = 0.17076884, "Z1" = 0.36590889, "eta~X1" = 1.1654944, "eta~X3" = 0.11297345, "Y2~eta" = 0.91569218, "Y2~X2" = 0.52324123, "Y3~eta" = 1.77781303, "Z1~eta" = 0.10836302, "Y1~~Y1" = 1.15857402, "eta~~eta" = 0.60432377, "Y2~~Y2" = 1.50666181, "Y3~~Y3" = 0.3155824, "Z1~~Z1" = 1.97827662, "Y1~~Z1" = 0.24216573, "Y2~~Z1" = 0.34997935)
    expect_equal(coef(eSSC1.lvm),GS, tol = 1e-6)
})

test_that("factor model - Cox correction", {
    eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")
    ## coef(eSSC2.lvm) - coef(e.lvm)
    
    GS <- c("eta" = 0.23991104, "Y2" = 0.29144175, "Y3" = 0.17832511, "Z1" = 0.36554264, "eta~X1" = 1.16545963, "eta~X3" = 0.11297009, "Y2~eta" = 0.90409861, "Y2~X2" = 0.52324125, "Y3~eta" = 1.75203851, "Z1~eta" = 0.10961229, "Y1~~Y1" = 1.20351764, "eta~~eta" = 0.5603901, "Y2~~Y2" = 1.47834635, "Y3~~Y3" = 0.30540849, "Z1~~Z1" = 1.99182361, "Y1~~Z1" = 0.25025483, "Y2~~Z1" = 0.3555143)
    expect_equal(coef(eSSC2.lvm),GS, tol = 1e-6)
})

## ** two factors model (correlation)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,
           eta1~eta2))

e.lvm <- estimate(m, d)

test_that("two factors model (correlation) - residuals correction", {
    eSSC1.lvm <- sCorrect(e.lvm, ssc = "residuals")

    GS <- c("eta1" = 0.1478569, "Y2" = 0.19515562, "Y3" = 0.37384111, "eta2" = 0.39767751, "Z2" = -0.19934296, "Z3" = -0.82545231, "eta1~eta2" = 0.36540063, "Y2~eta1" = 0.85064787, "Y3~eta1" = 0.88015853, "Y3~X1" = 1.29882077, "Z2~eta2" = 0.92947602, "Z3~eta2" = 1.95754764, "Z3~X3" = 1.22777389, "Y1~~Y1" = 0.82498637, "eta1~~eta1" = 2.19604428, "Y2~~Y2" = 1.67543559, "Y3~~Y3" = 0.90484948, "Z1~~Z1" = 1.20086385, "eta2~~eta2" = 0.7944319, "Z2~~Z2" = 1.41920933, "Z3~~Z3" = 0.21553652)
    expect_equal(coef(eSSC1.lvm),GS, tol = 1e-6)
}

test_that("two factors model (correlation) - Cox correction", {
    eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")

    GS <- c("eta1" = 0.15334368, "Y2" = 0.19799503, "Y3" = 0.37775788, "eta2" = 0.39767751, "Z2" = -0.18562638, "Z3" = -0.75078602, "eta1~eta2" = 0.35160357, "Y2~eta1" = 0.8409626, "Y3~eta1" = 0.86679841, "Y3~X1" = 1.29882077, "Z2~eta2" = 0.89498429, "Z3~eta2" = 1.76979176, "Z3~X3" = 1.22777389, "Y1~~Y1" = 0.88673054, "eta1~~eta1" = 2.22027337, "Y2~~Y2" = 1.7069357, "Y3~~Y3" = 0.87527537, "Z1~~Z1" = 1.26455176, "eta2~~eta2" = 0.74360185, "Z2~~Z2" = 1.43772669, "Z3~~Z3" = 0.31131063)

    expect_equal(coef(eSSC2.lvm), GS, tol = 1e-6)
}

## ** two factors model (covariance)
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1,eta1~X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3,eta2~X2,
           eta1~~eta2))

e.lvm <- estimate(m, d)

test_that("two factors model (correlation) - residuals correction", {
    eSSC1.lvm <- sCorrect(e.lvm, ssc = "residuals")

    GS <- c("eta1" = 0.1478569, "Y2" = 0.19515562, "Y3" = 0.37384111, "eta2" = 0.39767751, "Z2" = -0.19934296, "Z3" = -0.82545231, "eta1~eta2" = 0.36540063, "Y2~eta1" = 0.85064787, "Y3~eta1" = 0.88015853, "Y3~X1" = 1.29882077, "Z2~eta2" = 0.92947602, "Z3~eta2" = 1.95754764, "Z3~X3" = 1.22777389, "Y1~~Y1" = 0.82498637, "eta1~~eta1" = 2.19604428, "Y2~~Y2" = 1.67543559, "Y3~~Y3" = 0.90484948, "Z1~~Z1" = 1.20086385, "eta2~~eta2" = 0.7944319, "Z2~~Z2" = 1.41920933, "Z3~~Z3" = 0.21553652)
    expect_equal(coef(eSSC1.lvm),GS, tol = 1e-6)
}

test_that("two factors model (correlation) - Cox correction", {
    eSSC2.lvm <- sCorrect(e.lvm, ssc = "Cox")

    GS <- c("eta1" = 0.24295941, "Y2" = 0.18388462, "Y3" = 0.32464952, "eta2" = 0.37537296, "Z2" = -0.1851004, "Z3" = -0.77200208, "eta1~X1" = 1.08475547, "Y2~eta1" = 0.88925211, "Y3~eta1" = 1.12263979, "Y3~X1" = 0.82749017, "eta2~X2" = -0.10473803, "Z2~eta2" = 0.89399459, "Z3~eta2" = 1.82282964, "Z3~X3" = 1.21228718, "Y1~~Y1" = 0.98937181, "eta1~~eta1" = 0.9128264, "Y2~~Y2" = 1.59949563, "Y3~~Y3" = 0.75743488, "Z1~~Z1" = 1.29260368, "eta2~~eta2" = 0.72179897, "Z2~~Z2" = 1.46486865, "Z3~~Z3" = 0.22243345, "eta1~~eta2" = 0.14547098)

    expect_equal(coef(eSSC2.lvm), GS, tol = 1e-6)
    ## coef(e.lvm) - coef(eSSC2.lvm)
}

##----------------------------------------------------------------------
### test1-sCorrect-ssc.R ends here

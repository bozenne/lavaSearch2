### test-dVcov2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2018 (15:17) 
## Version: 
## Last-Updated: jan  3 2018 (18:41) 
##           By: Brice Ozenne
##     Update #: 5
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

context("dVcov2")
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

## ** lm
e.lvm <- estimate(lvm(Y1~X1+X2), data = dW)
e.lvm$prepareScore2 <- prepareScore2(e.lvm, second.order = TRUE, update = FALSE)
e.gls <- gls(Y1~X1+X2, data = dW, method = "ML")

dfVariance(e.lvm)
anova(e.gls)


test_that("linear regression: dVcov2",{
    ## lvm
    GS.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.lvm <- dVcov2(e.lvm, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.lvm, res.lvm)

    ## gls
    GS.gls <- dVcov2(e.lvm, adjust.residuals = FALSE,
                     numericDerivative = TRUE)
    res.gls <- dVcov2(e.lvm, adjust.residuals = FALSE,
                      numericDerivative = FALSE)
    expect_equal(GS.gls, res.gls)
    
})

test_that("linear regression: df adjusted",{
    df.adj.lvm <- dfVariance(e.lvm, adjust.residuals = TRUE)
    df.adj.gls <- dfVariance(e.gls, cluster = 1:n, adjust.residuals = TRUE,
                             numericDerivative = FALSE)
    GS <- c(rep(NROW(dW)-(n.param-1),n.param-1), (NROW(dW)-(n.param-1))/4)
    expect_equal(as.double(df.adj.gls),GS)
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


## *** lava - ok
expect_equal(as.double(logLik(e.lmer)),as.double(logLik(e.lvm)))

test_that("mixed model: df",{
    GS <- summary(e.lmer, ddf = "Satterthwaite")$coef[,"df"]
    
    df.lvm <- dfVariance(e.lvm, adjust.residuals = FALSE,
                         numericDerivative = TRUE)
    expect_equal(as.double(GS),
                 as.double(df.lvm[1:5]))

    df1.lme <- dfVariance(e.lme, adjust.residuals = FALSE,
                          numericDerivative = TRUE)
    expect_equal(GS, df1.lme[names(GS)])

    df2.lme <- dfVariance(e.lme, adjust.residuals = FALSE,
                          numericDerivative = FALSE)
    expect_equal(GS, df2.lme[names(GS)], tol = 1e-5)

    df1.gls <- dfVariance(e.gls, adjust.residuals = FALSE,
                          numericDerivative = TRUE)
    expect_equal(GS, df1.gls[names(GS)], tol = 1e-5)

    df2.gls <- dfVariance(e.gls, adjust.residuals = FALSE,
                          numericDerivative = FALSE)
    expect_equal(GS, df2.gls[names(GS)], tol = 1e-5)
})

test_that("mixed model: df adjusted",{
    GS <- summary(e.lmer, ddf = "Kenward-Roger")$coef[,"df"]
    ## get_Lb_ddf(e.lmer, c(0,1,0,0,0))
    ## get_Lb_ddf(e.lmer, c(0,0,0,1,0))
    
    df.adj.lvm <- dfVariance(e.lvm, adjust.residuals = TRUE,
                             numericDerivative = TRUE)
    df.adj.lvm

    df.adj.lme <- dfVariance(e.lme, adjust.residuals = TRUE,
                             numericDerivative = FALSE)
    df.adj.lme

    df.adj.gls <- dfVariance(e.gls, adjust.residuals = TRUE,
                             numericDerivative = FALSE)
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
    df.adj.lvm <- dfVariance(e.lvm, adjust.residuals = FALSE,
                             numericDerivative = TRUE)
    df.adj.lvm

    ## df.adj.lme <- dfVariance(e.lme,
    ##                          robust = FALSE, adjust.residuals = FALSE)

    system.time(
        df1.gls <- dfVariance(e.gls, adjust.residuals = FALSE,
                                  numericDerivative = TRUE)
    )
    system.time(
        df2.gls <- dfVariance(e.gls, adjust.residuals = FALSE,
                                  numericDerivative = FALSE)
    )
    system.time(
        df2.adj.gls <- dfVariance(e.gls, adjust.residuals = TRUE,
                                  numericDerivative = FALSE)
    )
    df2.adj.gls
})


##----------------------------------------------------------------------
### test-dVcov2.R ends here

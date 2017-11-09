### test-iid2-nlme.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: nov  6 2017 (12:57) 
## Version: 
## last-updated: nov  9 2017 (18:19) 
##           By: Brice Ozenne
##     Update #: 4
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

context("iid2-nlme")

n <- 5e1


## * mixed model
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
latent(mSim) <- ~eta
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- as.data.table(sim(mSim,n,latent = FALSE))
setkey(dW, "Id")
dL <- melt(dW,id.vars = c("G","Id"), variable.name = "time")
setkey(dL, "Id")

keep.cols <- c("eta","Y2","Y3","eta~G")

## ** gls
e.gls <- gls(value ~ time + G,
             correlation = corCompSymm(form =~ 1| Id),
             weight = varIdent(form = ~ 1|time),
             data = dL, method = "ML")
factor <- (e.gls$dims$N - e.gls$dims$p)/(e.gls$dims$N - e.gls$dims$p * (e.gls$method == "REML"))
name.param <- names(coef(e.gls))

test_that("gls: HC0/HC1", {
    iid2HC0.gls <- iid2(e.gls, return.df = FALSE, adjust.residuals = FALSE)

    VsandwichHC0.gls <- crossprod(iid2HC0.gls)[name.param,name.param]
    GS <- vcovCR(e.gls, type = "CR0", cluster = dL$Id) * factor^2
    expect_equal(as.double(GS),as.double(VsandwichHC0.gls), tolerance = 1e-10)
    
    GS <- vcovCR(e.gls, type = "CR1", cluster = dL$Id) * factor^2
    VsandwichHC1.gls <- crossprod(iid2HC0.gls)*n/(n-1)
    expect_equal(as.double(GS),as.double(VsandwichHC1.gls), tolerance = 1e-10)
})

test_that("gls: HC3", {
    iid2HC3.gls <- iid2(e.gls, return.df = FALSE, adjust.residuals = TRUE, power = 1)

    GS <- vcovCR(e.gls, type = "CR3", cluster = dL$Id) * factor^2    
    VsandwichHC3.gls <- crossprod(iid2HC3.gls)
    expect_equal(as.double(GS),as.double(VsandwichHC3.gls), tolerance = 1e-10) 
})

test_that("gls: HC2", {
    iid2HC2.gls <- iid2(e.gls, return.df = FALSE, adjust.residuals = TRUE, power = 0.5)

    GS <- vcovCR(e.gls, type = "CR2", cluster = dL$Id) * factor^2
    VsandwichHC2.gls <- crossprod(iid2HC2.gls)
    expect_equal(as.double(GS),as.double(VsandwichHC2.gls), tolerance = 1e-10)
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

test_that("lme: HC0/HC1", {
    iid2HC0.lme <- iid2(e.lme, return.df = FALSE, adjust.residuals = FALSE)
    iid2HC0.lvm <- iid2(e.lvm, return.df = FALSE, adjust.residuals = FALSE, use.information = TRUE)

    expect_equal(unname(iid2HC0.lme),unname(iid2HC0.lvm[,keep.cols]), tol = 1e-6)

    VsandwichHC0.lme <- crossprod(iid2HC0.lme)
    VsandwichHC0.lvm <- crossprod(iid2HC0.lvm[,keep.cols])
    GS <- vcovCR(e.lme, type = "CR0", cluster = dL$Id)
    expect_equal(as.double(GS),as.double(VsandwichHC0.lme), tolerance = 1e-10)
    expect_equal(as.double(GS),as.double(VsandwichHC0.lvm), tolerance = 1e-7)

    GS <- vcovCR(e.lme, type = "CR1", cluster = dL$Id)
    VsandwichHC1.lme <- crossprod(iid2HC0.lme)*n/(n-1)
    VsandwichHC1.lvm <- crossprod(iid2HC0.lvm[,keep.cols])*n/(n-1)
    expect_equal(as.double(GS),as.double(VsandwichHC1.lme), tolerance = 1e-10)
    expect_equal(as.double(GS),as.double(VsandwichHC1.lvm), tolerance = 1e-7)
})

test_that("lme: HC3", {
    iid2HC3.lme <- iid2(e.lme, return.df = FALSE, adjust.residuals = TRUE, power = 1)
    iid2HC3.lvm <- iid2(e.lvm, return.df = FALSE, adjust.residuals = TRUE, use.information = TRUE, power = 1)

    expect_equal(unname(iid2HC3.lme),unname(iid2HC3.lvm[,keep.cols]), tol = 1e-6)

    VsandwichHC3.lme <- crossprod(iid2HC3.lme)
    VsandwichHC3.lvm <- crossprod(iid2HC3.lvm[,keep.cols])
    GS <- vcovCR(e.lme, type = "CR3", cluster = dL$Id)
    expect_equal(as.double(GS),as.double(VsandwichHC3.lme), tolerance = 1e-10)
    expect_equal(as.double(GS),as.double(VsandwichHC3.lvm), tolerance = 1e-7)
})

test_that("lme: HC2", {
    iid2HC2.lme <- iid2(e.lme, return.df = FALSE, adjust.residuals = TRUE, power = 0.5)
    iid2HC2.lvm <- iid2(e.lvm, return.df = FALSE, adjust.residuals = TRUE, use.information = TRUE, power = 0.5)

    expect_equal(unname(iid2HC2.lme),unname(iid2HC2.lvm[,keep.cols]), tol = 1e-6)

    VsandwichHC2.lme <- crossprod(iid2HC2.lme)
    VsandwichHC2.lvm <- crossprod(iid2HC2.lvm[,keep.cols])
    GS <- vcovCR(e.lme, type = "CR2", cluster = dL$Id)
    expect_equal(as.double(GS),as.double(VsandwichHC2.lme), tolerance = 1e-10)
    expect_equal(as.double(GS),as.double(VsandwichHC2.lvm), tolerance = 1e-7)
})

#----------------------------------------------------------------------
### test-iid2-nlme.R ends here

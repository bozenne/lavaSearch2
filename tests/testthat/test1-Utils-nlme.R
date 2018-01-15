### test-Utils-nlme.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 16 2017 (10:36) 
## Version: 
## Last-Updated: jan 15 2018 (18:57) 
##           By: Brice Ozenne
##     Update #: 28
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
rm(list = ls(all.names = TRUE))
toRM <- names(sessionInfo()$otherPkgs)
if(!is.null(toRM)){
    lapply(paste('package:',,sep=""),
           detach,
           character.only=TRUE,unload=TRUE)
}
if(TRUE){ ## already called in test-all.R
    library(testthat)
    library(lavaSearch2)
    library(data.table)
    library(lava)    
}

library(nlme)
lava.options(symbols = c("~","~~"))

context("Utils-nlme")
n <- 5e1

.coef2 <- lavaSearch2:::.coef2
.coef2.gls <- lavaSearch2:::.coef2.gls
.coef2.lme <- lavaSearch2:::.coef2.lme

.getGroups2 <- lavaSearch2:::.getGroups2
.getGroups2.gls <- lavaSearch2:::.getGroups2.gls
.getGroups2.lme <- lavaSearch2:::.getGroups2.lme

## * data
mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,Y4~1*eta,eta~G+Gender))
latent(mSim) <- ~eta
categorical(mSim, labels = c("M","F")) <- ~Gender
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- as.data.table(sim(mSim,n,latent = FALSE))
setkey(dW, "Id")
dL <- melt(dW,id.vars = c("G","Id","Gender"), variable.name = "time")
setkey(dL, "Id")

## * Compound symmetry
e.lme <- lme(value ~ time + G + Gender,
             random = ~ 1|Id,
             data = dL,
             method = "ML")
e.lme.bis <- lme(value ~ time + G + Gender,
                 random = ~ 1|Id,
                 correlation = corCompSymm(),
                 data = dL,
                 method = "ML")
e.gls <- gls(value ~ time + G + Gender,
             correlation = corCompSymm(form=~ 1|Id),
             data = dL, method = "ML")

vecCoef.lme <- .coef2(e.lme)
vecCoef.lme.bis <- .coef2(e.lme.bis)
vecCoef.gls <- .coef2(e.gls)

groups.lme <- .getGroups2(e.lme)
groups.lme.bis <- .getGroups2(e.lme.bis)
groups.gls <- .getGroups2(e.gls)

test_that("Compound symmetry", {
    lsVcov.gls <- .getVarCov2(e.gls,
                              param = vecCoef.gls,
                              attr.param = attributes(vecCoef.gls),
                              endogenous = groups.gls$endogenous,
                              name.endogenous = groups.gls$name.endogenous,
                              n.endogenous = groups.gls$n.endogenous,
                              cluster = groups.gls$cluster,
                              n.cluster = groups.gls$n.cluster)

    expect_equal(unclass(getVarCov(e.gls)),
                 unname(lsVcov.gls$Omega))

    lsVcov.lme <- .getVarCov2(e.lme,
                              param = vecCoef.lme,
                              attr.param = attributes(vecCoef.lme),
                              endogenous = groups.lme$endogenous,
                              name.endogenous = groups.lme$name.endogenous,
                              n.endogenous = groups.lme$n.endogenous,
                              cluster = groups.lme$cluster,
                              n.cluster = groups.lme$n.cluster)

    expect_equal(unname(getVarCov(e.lme, type = "marginal", individuals = 1)[[1]]),
                 unname(lsVcov.lme$Omega))

    lsVcov.lme.bis <- .getVarCov2(e.lme.bis,
                                  param = vecCoef.lme.bis,
                                  attr.param = attributes(vecCoef.lme.bis),
                                  endogenous = groups.lme.bis$endogenous,
                                  name.endogenous = groups.lme.bis$name.endogenous,
                                  n.endogenous = groups.lme.bis$n.endogenous,
                                  cluster = groups.lme.bis$cluster,
                                  n.cluster = groups.lme.bis$n.cluster)

    expect_equal(unname(getVarCov(e.lme.bis, type = "marginal", individuals = 1)[[1]]),
                 unname(lsVcov.lme.bis$Omega))
})

## * Unstructured 
e.lme <- lme(value ~ time + G + Gender,
             random = ~ 1|Id,
             correlation = corSymm(),
             data = dL,
             method = "ML")
e.gls <- gls(value ~ time + G + Gender,
             correlation = corSymm(form=~ 1|Id),
             data = dL, method = "ML")

vecCoef.lme <- .coef2(e.lme)
vecCoef.gls <- .coef2(e.gls)

groups.lme <- .getGroups2(e.lme)
groups.gls <- .getGroups2(e.gls)

test_that("Unstructured ", {
    lsVcov.gls <- .getVarCov2(e.gls,
                              param = vecCoef.gls,
                              attr.param = attributes(vecCoef.gls),
                              endogenous = groups.gls$endogenous,
                              name.endogenous = groups.gls$name.endogenous,
                              n.endogenous = groups.gls$n.endogenous,
                              cluster = groups.gls$cluster,
                              n.cluster = groups.gls$n.cluster)

    expect_equal(unclass(getVarCov(e.gls)),
                 unname(lsVcov.gls$Omega))

    lsVcov.lme <- .getVarCov2(e.lme,
                              param = vecCoef.lme,
                              attr.param = attributes(vecCoef.lme),
                              endogenous = groups.lme$endogenous,
                              name.endogenous = groups.lme$name.endogenous,
                              n.endogenous = groups.lme$n.endogenous,
                              cluster = groups.lme$cluster,
                              n.cluster = groups.lme$n.cluster)

    expect_equal(unname(getVarCov(e.lme, type = "marginal", individuals = 1)[[1]]),
                 unname(lsVcov.lme$Omega))
})

## * Unstructured with weight
e.lme <- lme(value ~ time + G + Gender,
             random = ~ 1|Id,
             correlation = corSymm(),
             weight = varIdent(form = ~ 1|time),
             data = dL,
             method = "ML")
e.gls <- gls(value ~ time + G + Gender,
             correlation = corSymm(form=~ 1|Id),
             weight = varIdent(form = ~ 1|time),
             data = dL, method = "ML")

vecCoef.lme <- .coef2(e.lme)
vecCoef.gls <- .coef2(e.gls)

groups.lme <- .getGroups2(e.lme)
groups.gls <- .getGroups2(e.gls)

test_that("Unstructured ", {
    lsVcov.gls <- .getVarCov2(e.gls,
                              param = vecCoef.gls,
                              attr.param = attributes(vecCoef.gls),
                              endogenous = groups.gls$endogenous,
                              name.endogenous = groups.gls$name.endogenous,
                              n.endogenous = groups.gls$n.endogenous,
                              cluster = groups.gls$cluster,
                              n.cluster = groups.gls$n.cluster)

    expect_equal(unclass(getVarCov(e.gls)),
                 unname(lsVcov.gls$Omega))

    lsVcov.lme <- .getVarCov2(e.lme,
                              param = vecCoef.lme,
                              attr.param = attributes(vecCoef.lme),
                              endogenous = groups.lme$endogenous,
                              name.endogenous = groups.lme$name.endogenous,
                              n.endogenous = groups.lme$n.endogenous,
                              cluster = groups.lme$cluster,
                              n.cluster = groups.lme$n.cluster)

    expect_equal(unname(getVarCov(e.lme, type = "marginal", individuals = 1)[[1]]),
                 unname(lsVcov.lme$Omega))
})

## * 2 random effect model (error)
e.lme <- lme(value ~ time + G + Gender,
             random=~1|Id/Gender,
             data = dL,
             method = "ML")

expect_error(.getGroups2(e.lme))

## e.lme <- lme(value ~ time + G + Gender,
##              random=~1|Id,
##              correlation=corCompSymm(form = ~1|Gender),
##              data = dL,
##              method = "ML")
## incompatible

##----------------------------------------------------------------------
### test-Utils-nlme.R ends here

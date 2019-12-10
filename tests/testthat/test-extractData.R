### test-extractData.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: dec  9 2019 (09:31) 
## Version: 
## Last-Updated: dec  9 2019 (09:33) 
##           By: Brice Ozenne
##     Update #: 1
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

library(nlme)
library(lme4)
library(lava)
library(riskRegression)
library(survival)

lava.options(symbols = c("~","~~"))
context("Extract data from models")

## * settings
set.seed(10)
n <- 100

## * linear regression
Y1 <- rnorm(n, mean = 0)
Y2 <- rnorm(n, mean = 0.3)
df.sim <- rbind(data.frame(Y=Y1,G=1,Id = 1:5),
           data.frame(Y=Y2,G=2,Id = 1:5)
           )
m.lm <- lm(Y ~ G, data = df.sim)
test_that("extractData (lm)", {
    expect_named(extractData(m.lm, design.matrix = TRUE),
                 expected = c("(Intercept)","G"))
})

## * mixed model
test_that("extractData (gls/lme/lmer)", {
  m.gls <- gls(Y ~ G, weights = varIdent(form = ~ 1|Id), data = df.sim)
  expect_named(extractData(m.gls, design.matrix = TRUE),
               expected = c("Y","G","Id"))
  m.lme <- lme(Y ~ G, random = ~ 1|Id, data = df.sim)
  expect_named(extractData(m.lme, design.matrix = TRUE),
               expected = c("Y","G","Id"))
  m.lmer <- lmer(Y ~ G + (1|Id), data = df.sim)
  expect_named(extractData(m.lmer, design.matrix = TRUE),
               expected = c("(Intercept)","G"))
})

## * lvm
test_that("extractData (lvm)", {
  e <- estimate(lvm(Y ~ G), data = df.sim)
  expect_named(extractData(e, design.matrix = TRUE),
               expected = c("Y","G"))
})


## * Cox model
dt.sim <- sampleData(n, outcome = "survival")
test_that("extractData (survival)", {
    # no strata
    m.cox <- coxph(Surv(time, event) ~ X6 + X7, data = dt.sim, x = TRUE, y = TRUE)
    test.data <- extractData(m.cox, design.matrix = TRUE)
    expect_named(test.data,
                 expected = c("start","stop","status","X6","X7","strata"))
    # strata
    m.cox <- coxph(Surv(time, event) ~ strata(X1) + X6, data = dt.sim, x = TRUE, y = TRUE)
    expect_named(extractData(m.cox, design.matrix = TRUE),
                 expected = c("start","stop","status","X6","strata","X1"))
    expect_named(extractData(m.cox, design.matrix = FALSE),
                 expected = names(dt.sim))
})

######################################################################
### test-extractData.R ends here

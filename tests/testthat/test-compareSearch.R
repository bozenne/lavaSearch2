### test-compareSearch.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt  5 2017 (09:25) 
## Version: 
## last-updated: okt  5 2017 (12:05) 
##           By: Brice Ozenne
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(testthat)
library(data.table)
library(mvtnorm)

context("compareSearch")

## * example
n <- 100
m.sim <- lvm(Y~E+X1+0.5*X2+0.1*X3,~Z1+Z2+Z3)
m.base <- lvm(Y~E,~X1+X2+X3+Z1+Z2+Z3)

set.seed(10)
dt.sim <- as.data.table(sim(m.sim, n=100, latent = FALSE))
e.base <- estimate(m.base, data = dt.sim)

resCompare <- compareSearch(e.base, 
                            method.iid = "iid",
                            method.p.adjust = c("none","bonferroni","fdr","max"),
                            statistic = c("score","Wald"),trace = 5)
##----------------------------------------------------------------------
### test-compareSearch.R ends here

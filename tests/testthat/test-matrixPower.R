### test-matrixPower.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 24 2017 (09:21) 
## Version: 
## last-updated: jan 15 2018 (18:57) 
##           By: Brice Ozenne
##     Update #: 13
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
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

lava.options(symbols = c("~","~~"))

context("matrixPower")

## * tests
M <- matrix(rnorm(2e2),20,10)
Sigma <- var(M)

test_that("square root", {
    Sigma.half <- matrixPower(Sigma, power = 1/2, symmetric = FALSE)
    expect_true(all(abs(crossprod(Sigma.half)-Sigma) < 1e-12))
})
 
test_that("inverse", {
    Sigma.m1 <- matrixPower(Sigma, power = -1, symmetric = FALSE)
    expect_equal(Sigma.m1 %*% Sigma,diag(1,NROW(Sigma),NCOL(Sigma)))
})

#----------------------------------------------------------------------
### test-matrixPower.R ends here

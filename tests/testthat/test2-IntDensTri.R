### test-IntDensTri.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 31 2017 (16:32) 
## Version: 
## last-updated: jan 15 2018 (18:56) 
##           By: Brice Ozenne
##     Update #: 12
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

library(mvtnorm)
lava.options(symbols = c("~","~~"))

context("IntDensTri")

## * tests
# {{{ around 0
test_that("Integrate standard gaussian density (2D)", {
    p <- 2
    Sigma <- diag(p)
    mu <- rep(0, p)

    for(n in c(5,10,20,50,100)){ # n <- 5
        res <- IntDensTri(mu = mu, Sigma = Sigma, n=n, x.min=0)
        expect_equal(res$value,
                     1/2,
                     tol = 1e-6)
    }
})

test_that("Integrate standard gaussian density (3D)", {
    p <- 3
    Sigma <- diag(p)
    mu <- rep(0, p)

    for(n in c(5,10,20,50,100)){
        res <- IntDensTri(mu = mu, Sigma = Sigma, n=n, x.min=0, z.max=10)
        expect_equal(res$value,
                     1/2,
                     tol = 1e-6)
    }
})
# }}}

# {{{ far from 0
test_that("Integrate standard gaussian density (2D)", {
    p <- 2
    Sigma <- diag(p)
    mu <- c(10,0)

    for(n in c(5,10,20,50,100)){
        res <- IntDensTri(mu = mu, Sigma = Sigma, n=n, x.min=0)
        expect_equal(res$value,
                     1,
                     tol = 1e-6)
    }
})

test_that("Integrate standard gaussian density (3D)", {
    p <- 3
    Sigma <- diag(p)
    mu <- c(10,0,0)

    for(n in c(5,10,20,50,100)){
        res <- IntDensTri(mu = mu, Sigma = Sigma, n=n, x.min=0, z.max=10)
        expect_equal(res$value,
                     1,
                     tol = 1e-6)
    }
})
# }}}


#----------------------------------------------------------------------
### test-IntDensTri.R ends here

### test-IntDensTri.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 31 2017 (16:32) 
## Version: 
## last-updated: aug 31 2017 (16:40) 
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
library(mvtnorm)

context("IntDensTri")
# {{{ around 0
test_that("Integrate standard gaussian density (2D)", {
    p <- 2
    Sigma <- diag(p)
    mu <- rep(0, p)

    for(n in c(5,10,20,50,100)){
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

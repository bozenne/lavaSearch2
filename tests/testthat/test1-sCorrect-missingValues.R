### test1-sCorrect-missingValues.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  7 2018 (13:39) 
## Version: 
## Last-Updated: mar 13 2018 (13:25) 
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

## * header
rm(list = ls())
if(FALSE){ ## already called in test-all.R
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))

context("sCorrect: dealing with missing values")

## ** t test
## formula:
## df = \frac{ 2 * s_pool^2 }{ var(s_pool^2) }
##    = \frac{ ( s_X^2/m + s_Y^2/n )^2}{( s_X^4/(m(m-1)) + s_Y^4/(n(n-1)))}

## ## using the t test function
## e.ttest <- t.test(dW$Y1,dW$Y2)
## e.ttest$parameter

## ## by hand
## sX1 <- var(dW$Y1)/n
## sX2 <- var(dW$Y2)/n
## df <- (sX1+sX2)^2/(sX1^2/(n-1) + sX2^2/(n-1))

## df-e.ttest$parameter


##----------------------------------------------------------------------
### test1-sCorrect-missingValues.R ends here

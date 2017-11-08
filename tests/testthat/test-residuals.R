### test-residuals.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov  8 2017 (09:08) 
## Version: 
## Last-Updated: nov  8 2017 (15:04) 
##           By: Brice Ozenne
##     Update #: 6
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

context("residuals2")

n <- 5e1

mSim <- lvm(c(Y1~1*eta,Y2~1*eta,Y3~1*eta,eta~G))
latent(mSim) <- ~eta
transform(mSim,Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
dW <- as.data.table(sim(mSim,n,latent = FALSE))
setkey(dW, "Id")
dL <- melt(dW,id.vars = c("G","Id"), variable.name = "time")
setkey(dL, "Id")
dL$Z1 <- rnorm(NROW(dL))

## * several linear models
m <- lvm(Y~X+G,G~X)

## * mixed model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~1*eta1,eta1~beta*G1,
           Z1~eta2,Z2~eta2,Z3~1*eta2,eta2~beta*G2)
         )
latent(m) <- ~eta1+eta2
e.lvm <- estimate(m,sim(m, 5e1))
e.lvm2 <- e.lvm
e.lvm2$prepareScore2 <- prepareScore2(e.lvm2)

test_that("equivalence residuals2.lvm residuals.lvm", {
    test <- residuals2(e.lvm, adjust.residuals = FALSE)
    test2 <- residuals2(e.lvm2, adjust.residuals = FALSE)    
    GS <- residuals(e.lvm)
    expect_equal(GS,test)
    expect_equal(GS,test2)
})



##----------------------------------------------------------------------
### test-residuals.R ends here

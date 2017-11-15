### test-iid2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (13:31) 
## Version: 
## last-updated: nov 15 2017 (13:45) 
##           By: Brice Ozenne
##     Update #: 116
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

context("iid2-lava")

n <- 5e1

## * linear model
p <- 3
X.name <- paste0("X",1:p)
link.lvm <- paste0("Y~",X.name)
formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))

set.seed(10)
m <- lvm(formula.lvm)
transform(m,Id~Y) <- function(x){1:NROW(x)}
set.seed(10)
d <- sim(m,n)

e.lm <- lm(formula.lvm,data=d)
e.lvm <- estimate(lvm(formula.lvm),data=d)

## ** iid2 matches iid
test_that("iid2 matches iid", {
    e.iid2.lm <- iid2(e.lm, adjust.residuals = FALSE)
    GS1 <- iid(e.lm)
    attr(GS1, "bread") <- NULL
    expect_equal(e.iid2.lm, GS1)

    e1.iid2.lvm <- iid2(e.lvm, adjust.residuals = FALSE)
    expect_equal(unname(e1.iid2.lvm[,1:4]), unname(GS1))
})


## ** iid2 lvm matches iid2 lm
test_that("iid2 lvm matches iid2 lm", {
    for(iAdj in c(FALSE,TRUE)){ # iAdj <- 1
        for(iPower in c(0.5,1)){ # iPower <- 1
            e.iid2.lm <- iid2(e.lm, adjust.residuals = iAdj, power = iPower)
            e0.iid2.lvm <- iid2(e.lvm, adjust.residuals = iAdj, power = iPower, use.information = TRUE)
            expect_equal(unname(e.iid2.lm), unname(e0.iid2.lvm[,1:4]), tolerance = 1e-10)
        }
    }
})

## ** iid2 matches clubSandwich
test_that("iid2.lm matches clubSandwich", {
    eHC2.iid2 <- iid2(e.lm, adjust.residuals = TRUE, power = 0.5)
    VsandwichHC2.lm <- crossprod(eHC2.iid2)
    expect_equal(as.double(vcovCR(e.lm, type = "CR2", cluster = d$Id)),
                 as.double(VsandwichHC2.lm))

    eHC3.iid2 <- iid2(e.lm, adjust.residuals = TRUE, power = 1)
    VsandwichHC3.lm <- crossprod(eHC3.iid2)
    expect_equal(as.double(vcovCR(e.lm, type = "CR3", cluster = d$Id)),
                 as.double(VsandwichHC3.lm))
})

test_that("iid2.lvm matches clubSandwich", {
    eHC2.iid2 <- iid2(e.lvm, adjust.residuals = TRUE, power = 0.5)
    VsandwichHC2.lvm <- crossprod(eHC2.iid2)[1:4,1:4]
    expect_equal(as.double(vcovCR(e.lm, type = "CR2", cluster = d$Id)),
                 as.double(VsandwichHC2.lvm))

    eHC3.iid2 <- iid2(e.lvm, adjust.residuals = TRUE, power = 1)
    VsandwichHC3.lvm <- crossprod(eHC3.iid2)[1:4,1:4]
    expect_equal(as.double(vcovCR(e.lm, type = "CR3", cluster = d$Id)),
                 as.double(VsandwichHC3.lvm))
})




#----------------------------------------------------------------------
### test-iid2.R ends here

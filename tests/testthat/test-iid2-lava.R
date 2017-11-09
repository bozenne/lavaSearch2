### test-iid2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (13:31) 
## Version: 
## last-updated: nov  9 2017 (18:19) 
##           By: Brice Ozenne
##     Update #: 108
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
    e.iid2.lm <- iid2(e.lm, return.df = FALSE, adjust.residuals = FALSE)
    GS1 <- iid(e.lm)
    attr(GS1, "bread") <- NULL
    expect_equal(e.iid2.lm, GS1)

    e1.iid2.lvm <- iid2(e.lvm, return.df = FALSE, adjust.residuals = FALSE)
    expect_equal(unname(e1.iid2.lvm[,1:4]), unname(GS1))
})


tempo <- iid2(e.lvm, return.df = TRUE, adjust.residuals = FALSE)
iid2(e.lvm, p = pars(e.lvm))[1,]
score2(e.lvm, p = pars(e.lvm) + c(0,0,0,0,1))[1,]

solve(vcov(e.lvm))
n*(1/coef(e.lvm)["Y~~Y"])^2/2

crossprod(model.matrix(e.lm))/coef(e.lvm)["Y~~Y"]

## ** iid2 lvm matches iid2 lm
test_that("iid2 lvm matches iid2 lm", {
    for(iAdj in c(FALSE,TRUE)){ # iAdj <- 1
        for(iPower in c(0.5,1)){ # iPower <- 1
            e.iid2.lm <- iid2(e.lm, return.df = FALSE, adjust.residuals = iAdj, power = iPower)
            e0.iid2.lvm <- iid2(e.lvm, return.df = FALSE, adjust.residuals = iAdj, power = iPower, use.information = TRUE)
            expect_equal(unname(e.iid2.lm), unname(e0.iid2.lvm[,1:4]), tolerance = 1e-10)
        }
    }
})

## ** iid2 matches clubSandwich
test_that("iid2.lm matches clubSandwich", {
    eHC2.iid2 <- iid2(e.lm, return.df = FALSE, adjust.residuals = TRUE, power = 0.5)
    VsandwichHC2.lm <- crossprod(eHC2.iid2)
    expect_equal(as.double(vcovCR(e.lm, type = "CR2", cluster = d$Id)),
                 as.double(VsandwichHC2.lm))

    eHC3.iid2 <- iid2(e.lm, return.df = FALSE, adjust.residuals = TRUE, power = 1)
    VsandwichHC3.lm <- crossprod(eHC3.iid2)
    expect_equal(as.double(vcovCR(e.lm, type = "CR3", cluster = d$Id)),
                 as.double(VsandwichHC3.lm))
})

test_that("iid2.lvm matches clubSandwich", {
    eHC2.iid2 <- iid2(e.lvm, return.df = FALSE, adjust.residuals = TRUE, power = 0.5)
    VsandwichHC2.lvm <- crossprod(eHC2.iid2)[1:4,1:4]
    expect_equal(as.double(vcovCR(e.lm, type = "CR2", cluster = d$Id)),
                 as.double(VsandwichHC2.lvm))

    eHC3.iid2 <- iid2(e.lvm, return.df = FALSE, adjust.residuals = TRUE, power = 1)
    VsandwichHC3.lvm <- crossprod(eHC3.iid2)[1:4,1:4]
    expect_equal(as.double(vcovCR(e.lm, type = "CR3", cluster = d$Id)),
                 as.double(VsandwichHC3.lvm))
})


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


#----------------------------------------------------------------------
### test-iid2.R ends here

### test-iid2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (13:31) 
## Version: 
## last-updated: nov  7 2017 (19:31) 
##           By: Brice Ozenne
##     Update #: 106
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
# e.iid2.lm <- iid2(e.lvm, use.information = TRUE)
test_that("iid2 matches iid", {
    e.iid2.lm <- iid2(e.lm, return.df = FALSE, adjust.residuals = FALSE)
    GS1 <- iid(e.lm)
    attr(GS1, "bread") <- NULL
    expect_equal(e.iid2.lm, GS1)

    e1.iid2.lvm <- iid2(e.lvm, return.df = FALSE, adjust.residuals = FALSE)
    expect_equal(unname(e1.iid2.lvm[,1:4]), unname(GS1))
})

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

## ** gls
e.gls <- gls(value ~ time + G,
             correlation = corCompSymm(form =~ 1| Id),
             weight = varIdent(form = ~ 1|time),
             data = dL, method = "ML")
factor <- (e.gls$dims$N - e.gls$dims$p)/(e.gls$dims$N - e.gls$dims$p * (e.gls$method == "REML"))

test_that("gls: HC0/HC1", {
    iid2HC0.gls <- iid2(e.gls, return.df = FALSE, adjust.residuals = FALSE)

    VsandwichHC0.gls <- crossprod(iid2HC0.gls)
    GS <- vcovCR(e.gls, type = "CR0", cluster = dL$Id) * factor^2
    expect_equal(as.double(GS),as.double(VsandwichHC0.gls), tolerance = 1e-10)
    
    GS <- vcovCR(e.gls, type = "CR1", cluster = dL$Id) * factor^2
    VsandwichHC1.gls <- crossprod(iid2HC0.gls)*n/(n-1)
    expect_equal(as.double(GS),as.double(VsandwichHC1.gls), tolerance = 1e-10)
})

test_that("gls: HC3", {
    iid2HC3.gls <- iid2(e.gls, return.df = FALSE, adjust.residuals = TRUE, power = 1)

    GS <- vcovCR(e.gls, type = "CR3", cluster = dL$Id) * factor^2    
    VsandwichHC3.gls <- crossprod(iid2HC3.gls)
    expect_equal(as.double(GS),as.double(VsandwichHC3.gls), tolerance = 1e-10) 
})

test_that("gls: HC2", {
    iid2HC2.gls <- iid2(e.gls, return.df = FALSE, adjust.residuals = TRUE, power = 0.5)

    GS <- vcovCR(e.gls, type = "CR2", cluster = dL$Id) * factor^2
    VsandwichHC2.gls <- crossprod(iid2HC2.gls)
    expect_equal(as.double(GS),as.double(VsandwichHC2.gls), tolerance = 1e-10)
})


## ** lme

m <- lvm(c(Y~X+G,G~X))
dW <- sim(m, 1e1)
e.lvm <- estimate(m, dW)

score2(e.lvm, return.df = FALSE, adjust.residuals = FALSE, return.vcov.param = TRUE)
head(score(e.lvm, indiv = TRUE))

vcov(e.lvm)
prepareScore2(e.lvm)$dmu.dtheta[6]

m00 <- lm(Y~X+G, data = dW)
m0 <- estimate(lvm(Y~X+G), data = dW)
score(m0, indiv = TRUE)
score2(m0, indiv = TRUE, adjust.residuals = FALSE)

dW$Y-coef(m0)[1]-coef(m0)[2]*dW$X-coef(m0)[3]*dW$G
dW$Y-predict(m0)




residuals(m0)
residuals(e.lvm)
predict(e.lvm)
predictlvm(e.lvm)

residuals(e.lvm, data = )
lava:::predict.lvmfit
dW$G
dW$X
coefType(m)

m <- lvm(c(Y1~eta,Y2~eta,Y3~eta+X1))
covariance(m) <- Y1~Y2
dW <- sim(m, 1e2)
e.lvm <- estimate(m, dW)

## ** iid2 matches iid
test_that("iid2 matches iid", {
    e1.iid2.lvm <- iid2(e.lvm, return.df = FALSE, adjust.residuals = FALSE)
    e1.iid2.lvm <- score2(e.lvm, return.df = FALSE, adjust.residuals = FALSE)

    
    GS1 <- iid(e.lvm)
    attr(GS1, "bread") <- NULL

    expect_equal(unname(e1.iid2.lvm[,1:4]), unname(GS1))
    
})



#----------------------------------------------------------------------
### test-iid2.R ends here

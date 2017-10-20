### test-iid2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (13:31) 
## Version: 
## last-updated: okt 19 2017 (19:04) 
##           By: Brice Ozenne
##     Update #: 30
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(testthat)
library(sandwich)

context("iid2")

n <- 5e1

## * linear model
p <- 3
X.name <- paste0("X",1:p)
link.lvm <- paste0("Y~",X.name)
formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))

set.seed(10)
m <- lvm(formula.lvm)
distribution(m,~Id) <- sequence.lvm(0)
set.seed(10)
d <- sim(m,n)

e.lm <- lm(formula.lvm,data=d)
e.lvm <- estimate(lvm(formula.lvm),data=d)

## ** iid2 matches iid
test_that("iid2 matches iid", {
    e.iid2.lm <- iid2(e.lm, return.df = FALSE, adjust.residuals = FALSE, use.information = FALSE, Dmethod = lava.options()$Dmethod)
    GS1 <- iid(e.lm)
    attr(GS1, "bread") <- NULL
    expect_equal(e.iid2.lm, GS1)

    e.iid2.lvm <- iid2(e.lvm, return.df = FALSE, adjust.residuals = FALSE, use.information = FALSE, Dmethod = lava.options()$Dmethod)
    GS2 <- iid(e.lvm)
    attr(GS2, "bread") <- NULL
    expect_equal(e.iid2.lvm, GS2)
})

## ** iid2 lvm matches iid2 lm
test_that("iid2 lvm matches iid2 lm", {
    for(iAdj in c(FALSE,TRUE)){ # iAdj <- 1
    e.iid2.lm <- iid2(e.lm, return.df = FALSE, adjust.residuals = iAdj)
    
    e0.iid2.lvm <- iid2(e.lvm, use.information = TRUE, adjust.residuals = iAdj)
    expect_equal(unname(e.iid2.lm), unname(e0.iid2.lvm[,1:4]), tolerance = 1e-10)

    iData <- model.matrix(e.lm)
    XX_m1 <- solve(t(data)%*%data)
    H <- data %*% XX_m1 %*% t(data)    

    iEpsilon <- residuals(e.lm)/(1-diag(iData)^(1/2))
    sweep(iData, MARGIN = 1, FUN = "*", STATS = iEpsilon)/score(e.lvm, indiv = TRUE)[,1:4]
    sweep(iData, MARGIN = 1, FUN = "*", STATS = iEpsilon)/score2(e.lvm, indiv = TRUE)[,1:4]
    
    
    e1.iid2.lvm <- iid2(e.lvm, use.information = FALSE, adjust.residuals = iAdj, Dmethod = "simple")
    expect_equal(e0.iid2.lvm, e1.iid2.lvm, tolerance = 1e-4)
    expect_equal(unname(e.iid2.lm), unname(e1.iid2.lvm[,1:4]), tolerance = 1e-4)
    e2.iid2.lvm <- iid2(e.lvm, use.information = FALSE, adjust.residuals = iAdj, Dmethod = "Richardson")
    expect_equal(e0.iid2.lvm, e2.iid2.lvm, tolerance = 1e-10)
    expect_equal(unname(e.iid2.lm), unname(e2.iid2.lvm[,1:4]), tolerance = 1e-10)
    }
}

test_that("iid2 matches sandwich", {
    eHC2.iid2 <- iid2(e.lm, return.df = FALSE, adjust.residuals = TRUE, alpha = 0.5)
    VsandwichHC2.lm <- t(eHC2.iid2) %*% eHC2.iid2
    expect_equal(unname(vcovCL(e.lm, type = "HC2")),unname(VsandwichHC2.lm))

    eHC3.iid2 <- iid2(e.lm, return.df = FALSE, adjust.residuals = TRUE, alpha = 1)
    VsandwichHC3.lm <- t(eHC3.iid2) %*% eHC3.iid2
    expect_equal(unname(vcovCL(e.lm, type = "HC3")),unname(VsandwichHC3.lm))

    eHC3.iid2.lvm <- iid2(e.lvm, alpha = 1, use.information = TRUE, adjust.residuals = TRUE)
    VsandwichHC3.lvm <- crossprod(eHC3.iid2.lvm)
    expect_equal(unname(VsandwichHC3.lvm[1:4,1:4]), unname(VsandwichHC3.lm))
    expect_equal(unname(vcovCL(e.lm, type = "HC3")),unname(Vsandwich0.lvm[1:4,1:4]))
})

score2(e.lvm, indiv = TRUE)
debug(score)
score(e.lvm, indiv = TRUE)

## ** simulation
if(FALSE){
    warper <- function(n){

        d <- data.table(Y = rnorm(n),X = rnorm(n), Id = 1:n)

        e.lm <- lm(Y~X,data=d)
        e.iid <- sqrt(n) * iid(e.lm)
        e.iid2 <- sqrt(n) * iid2(e.lm)

        sd.lm <- sqrt(vcov(e.lm)[2,2])
        sdR.iid <- sd(e.iid[,2])
        sdR.iid2 <- sd(e.iid2[,2])

        Wald.lm <- coef(e.lm)[2]/sd.lm
        Wald.iid <- coef(e.lm)[2]/sdR.iid
        Wald.iid2 <- coef(e.lm)[2]/sdR.iid2

    
        out <- NULL
    
        out <- rbind(out, data.table(p = summary(e.lm)$coef[2,4], method = "lm"))
        out <- rbind(out, data.table(p = 2*(1-pt(abs(Wald.lm), df = n-2)), method = "information"))
        out <- rbind(out, data.table(p = 2*(1-pt(abs(Wald.iid), df = n-2)), method = "iid"))
        out <- rbind(out, data.table(p = 2*(1-pt(abs(Wald.iid2), df = n-2)), method = "iid2"))
        out[, n:= n]
        return(out)
    
    }

    n.rep <- 1000
    resApply <- pbapply::pblapply(1:n.rep, function(x){
        dt <- warper(20)[, iRep := x]
        return(dt)
    })

    dt.pval <- rbindlist(resApply)
    dt.pval[,mean(p <= 0.05),by = "method"]
}


## * lvm

#----------------------------------------------------------------------
### test-iid2.R ends here

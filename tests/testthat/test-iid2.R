### test-iid2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (13:31) 
## Version: 
## last-updated: okt 12 2017 (15:24) 
##           By: Brice Ozenne
##     Update #: 18
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

m <- lvm(formula.lvm)
distribution(m,~Id) <- sequence.lvm(0)
set.seed(10)
d <- sim(m,n)

e.lm <- lm(formula.lvm,data=d)
e.lvm <- estimate(lvm(formula.lvm),data=d)


test_that("iid2 match sandwich", {
    e.iid2.lm <- iid2(e.lm)
    expect_equal(e.iid2.lm, iid2(e.lm, data = d))

    e.iid2.lvm <- iid2(e.lvm)
    
    Vsandwich.manual <- t(e.iid2.lm) %*% e.iid2.lm
    expect_equal(unname(vcovCL(e.lm, type = "HC3")),unname(Vsandwich.manual))
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

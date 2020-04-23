### test-mmm2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 29 2017 (15:22) 
## Version: 
## Last-Updated: apr 23 2020 (12:53) 
##           By: Brice Ozenne
##     Update #: 126
##----------------------------------------------------------------------
## 
### Commentary: 
## Test battery:
##  - list of linear regressions:
##      Compare multcomp:::glht to lavaSearch2:::glht2 (no small sample adjustment)
##      * standard error derived from the information matrix
##      * robust standard error derived from the iid decomposition
##
##  - latent variable models:
##      Compare multcomp:::glht to lavaSearch2:::glht2 (no small sample adjustment)
##      * standard error derived from the information matrix, with or without second member
##
##  - list of latent variable model (for linear regressions):
##      Compare multcomp:::glht to lavaSearch2:::glht2 (no small sample adjustment)
##      * standard error derived from the information matrix
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

library(multcomp)
library(sandwich)
library(emmeans)
lava.options(symbols = c("~","~~"))

context("multcomp - mmm")

## * simulation
mSim <- lvm(c(Y1,Y2,Y3,Y4)~ beta * eta,
            E ~ 1, Y1 ~ 0.25*T1 + 0.5*T2 + 0.05*T3)
latent(mSim) <- "eta"
set.seed(10)
n <- 1e2

df.data <- lava::sim(mSim, n, latent = FALSE, p = c(beta = 1))
df.data$eY1 <- exp(df.data$Y1)
## * linear regressions with logical constrains
e.lm <- lm(Y1 ~ T1 + T2 + T3, data = df.data)
## summary(e.lm)

test_that("glht vs. glht2 (logical constrains)", {
    e.glht <- glht(e.lm, linfct = c("T2-T1=0",
                                    "T2-T3=0",
                                    "T1-T3=0"))
    summary(e.glht, test = adjusted("none"))
    summary(e.glht, test = adjusted("bonferroni"))

    e.glht2 <- glht2(e.lm, linfct = c("T2-T1=0",
                                      "T2-T3=0",
                                      "T1-T3=0"))
    
    expect_equal(e.glht$vcov,e.glht2$vcov[rownames(e.glht$vcov),colnames(e.glht$vcov)], tol = 1e-6)
    expect_equal(e.glht$coef,e.glht2$coef[names(e.glht$coef)], tol = 1e-6)

    eS.glht <- summary(e.glht, test = adjusted("Shaffer"))
    eS.glht2 <- summary(e.glht2, test = adjusted("Shaffer"))

    expect_equivalent(eS.glht$test[c("coefficients","sigma","tstat","pvalues")],
                      eS.glht2$test[c("coefficients","sigma","tstat","pvalues")], tol = 1e-6)
})

test_that("glht2 (back-transformation)", {
    e.log.lm <- lm(log(eY1) ~ T1 + T2 + T3, data = df.data)

    e.glht2 <- glht2(e.log.lm, linfct = c("T1","T2","T3"))
    df.glht2 <- summary(e.glht2, transform = "exp", test = adjusted("none"))$table2

    e.glht2.bis <- glht2(e.log.lm, linfct = c("T3"))
    df.glht2.bis <- summary(e.glht2.bis, transform = "exp", test = adjusted("none"))$table2

    expect_equal(as.double(df.glht2[3,]) , as.double(df.glht2.bis[1,]))
})

## * list of linear regressions
name.Y <- setdiff(endogenous(mSim),"E")[1:2]
n.Y <- length(name.Y)

ls.lm <- lapply(name.Y, function(iY){
    eval(parse( text = paste("lm(",iY,"~E, data = df.data)")))
})
names(ls.lm) <- name.Y
class(ls.lm) <- "mmm"

test_that("glht vs. glht2 (list lm): information std", {
    e.glht <- glht(ls.lm, mlf("E = 0"))

    resC <- createContrast(ls.lm, linfct = "E", add.variance = TRUE)
    name.all <- colnames(resC$contrast)
    name.mean <- name.all[-grep("sigma",name.all)]

    e.glht2 <- glht2(ls.lm, linfct = resC$contrast,
                     ssc = "residuals", robust = FALSE)

    expect_equal(e.glht$vcov,
                 e.glht2$vcov[name.mean,name.mean], tol = 1e-6)
    expect_equal(e.glht$coef,e.glht2$coef[names(e.glht$coef)], tol = 1e-6)
    expect_equivalent(e.glht$linfct,e.glht2$linfct[,name.mean])

    eS.glht <- summary(e.glht)
    eS.glht2 <- summary(e.glht2)

    expect_equal(eS.glht$test$tstat, eS.glht2$test$tstat, tol = 1e-6)
})
     
test_that("glht vs. glht2 (list ml): robust std", {
    e.glht <- glht(ls.lm, mlf("E = 0"), vcov = sandwich)

    resC <- createContrast(ls.lm, linfct = "E", add.variance = TRUE)
    name.all <- colnames(resC$contrast)
    name.mean <- name.all[-grep("sigma",name.all)]
    e.glht2 <- glht2(ls.lm, linfct = resC$contrast,
                     ssc = NA, robust = TRUE)

    expect_equivalent(e.glht$vcov,
                      e.glht2$vcov[name.mean,name.mean])
    expect_equal(e.glht$coef,e.glht2$coef[name.mean])
    expect_equivalent(e.glht$linfct,e.glht2$linfct[,name.mean])

    eS.glht <- summary(e.glht)
    eS.glht2 <- summary(e.glht2)

    expect_equal(eS.glht$test$tstat, eS.glht2$test$tstat)
    ## cannot compare p.values
    ## because some are based on a student law and others on a gaussian law
})

test_that("glht2 vs. lava (ml): robust std", {
    lsRed.lm <- ls.lm[1:2]
    class(lsRed.lm) <- "mmm"

    resC <- createContrast(lsRed.lm, linfct = "E", add.variance = TRUE)
    name.all <- colnames(resC$contrast)
    name.mean <- name.all[-grep("sigma",name.all)]
    e.glht2 <- glht2(lsRed.lm, linfct = resC$contrast,
                     ssc = NA, robust = TRUE, df = NA)

    GS <- estimate(ls.lm[[1]], cluster = 1:n)$coefmat
    test <- summary(e.glht2, test = adjusted("none"))$test
    
    expect_equal(as.double(test$sigma[1]), GS["E","Std.Err"], tol = 1e-8)
    expect_equal(as.double(test$pvalues[1]), GS["E","P-value"], tol = 1e-8)
    ##
})

test_that("glht vs. calcDistMaxIntegral", {
    e.glht <- glht(ls.lm, mlf("E = 0"), vcov = sandwich)
    res.GS <- summary(e.glht)

    iid.tempo <- do.call(cbind,lapply(ls.lm, iid)) %*% t(e.glht$linfct)
    beta <- unlist(lapply(ls.lm, coef)) %*% t(e.glht$linfct)
    beta.var <- diag(crossprod(iid.tempo))
    z.value <- beta/sqrt(beta.var)
    res.Search <- calcDistMaxIntegral(as.vector(z.value),
                                      iid = iid.tempo, quantile.compute = FALSE,
                                      df = NULL, trace = FALSE, alpha = 0.05)
    
    expect_equal(as.double(res.Search$p.adjust),
                 as.double(res.GS$test$pvalues),
                 tolerance = attr(res.GS$test$pvalues, "error")
                 )
    ## cannot compare p.values
    ## because some are based on a student law and others on a gaussian law
})



## * lvm
## names(df.data)

m.lvm <- lvm(c(Y1,Y2,Y3)~ eta, eta ~ E,
             Y1~Y4)
e.lvm <- estimate(m.lvm, df.data)

test_that("glht vs. glht2 (lvm): information std", {

    resC <- createContrast(e.lvm, linfct = c("eta~E","Y2=1","Y3=1"))
    e.glht.null <- glht(e.lvm, linfct = resC$contrast)
    e.glht.H1 <- glht(e.lvm, linfct = resC$contrast, rhs = resC$null)

    e.glht2.null <- glht2(e.lvm, linfct = resC$contrast, rhs = rep(0,3),
                          ssc = NA)
    e.glht2.H1 <- glht2(e.lvm, linfct = resC$contrast, rhs = resC$null,
                        ssc = NA)


    eS.glht.null <- summary(e.glht.null)
    eS.glht.H1 <- summary(e.glht.H1)
    eS.glht2.null <- summary(e.glht2.null)
    eS.glht2.H1 <- summary(e.glht2.H1)

    expect_equal(unname(eS.glht.null$test$tstat),
                 unname(eS.glht2.null$test$tstat))
    expect_equal(unname(eS.glht.H1$test$tstat),
                 unname(eS.glht2.H1$test$tstat))
    ## cannot compare p.values
    ## because some are based on a student law and others on a gaussian law
})

## * list of lvm
mmm.lvm <- mmm(Y1 = estimate(lvm(Y1~E), data = df.data),
               Y2 = estimate(lvm(Y2~E), data = df.data),
               Y3 = estimate(lvm(Y3~E), data = df.data),
               Y4 = estimate(lvm(Y4~E), data = df.data)
               )

test_that("glht2 (list lvm): information std", {

    ##    
    resC <- createContrast(mmm.lvm, linfct = paste0("Y",1:4,": Y",1:4,"~E"))
    lvm2.glht <- glht2(mmm.lvm, linfct = resC$contrast,
                       ssc = NA, robust = FALSE)
    lvm2.sglht <- summary(lvm2.glht)    

    lvm3.glht <- glht2(mmm.lvm, linfct = resC$contrast,
                       rhs = rnorm(4),
                       ssc = NA, robust = FALSE)
    lvm3.sglht <- summary(lvm3.glht)    
})

##----------------------------------------------------------------------
### test-mmm2.R ends here

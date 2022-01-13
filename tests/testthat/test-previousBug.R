### test-previousBug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 19 2019 (10:17) 
## Version: 
## Last-Updated: Jan 12 2022 (14:55) 
##           By: Brice Ozenne
##     Update #: 29
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
## rm(list = ls(all.names = TRUE))
if(TRUE){ ## already called in test-all.R
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))

context("Previous bugs")


## * Brice, 01/23/20 1:41 , sCorrect
## keep.row <- sample.int(139,30)
## dd <- df.longiHCENRER[keep.row,c("neocortex.log","response.w4","age","sex","sert","sb.per.kg","group")]
## dd$sex <- as.numeric(dd$sex)
## dd$group <- as.numeric(dd$group)
## dd$sert <- as.numeric(dd$sert)
## dd$neocortex.log <- dd$neocortex.log+rnorm(NROW(dd), sd = 0.1)

ddW <- data.frame("neocortex.log" = c(-0.5114974, -0.8249681, -0.3910322, -0.2941140, -0.4850710, -0.7354951, -0.3005508, -0.3956807, -0.5279078, -0.4705333, -0.4072993, -0.5494491, -0.9464054, -0.2632149, -0.2781923, -1.0308939, -0.5637023, -0.2689296, -0.4375162, -0.5351115, -0.8410153, -0.3022184, -0.5711485, -0.5155027, -0.2262801, -0.4262566, -0.4353830, -0.7079919, -0.1529665, -0.3954827), 
                 "group" = c(rep(1,10),rep(2,20)),
                 "id" = 1:30)
ddW$group <- factor(ddW$group, levels = 1:2)

e.GS <- gls(neocortex.log ~ group-1, method = "REML", weight = varIdent(form =~1|group), data = ddW)
## same as
e.lm1 <- lm(neocortex.log ~ 1, ddW[ddW$group==1, ]) 
e.lm2 <- lm(neocortex.log ~ 1, ddW[ddW$group==2, ]) 

ddL <- reshape2::dcast(ddW, value.var = "neocortex.log", id~group)
names(ddL) <- c("id","G1","G2")
m <- lvm(G1~1,G2~1)
e.lvm <- estimate(m, ddL)

test_that("sCorrect in stratified GLS equivalent to separate LM", {
    eSSC.res <- estimate2(e.lvm)

    ## eSSC.cox <- estimate2(e.lvm, ssc = "Cox")

    GS <- c(mu1 = mean(ddL$G1, na.rm = TRUE),
            mu2 = mean(ddL$G2, na.rm = TRUE),
            sigma1 = var(ddL$G1, na.rm = TRUE),
            sigma2 = var(ddL$G2, na.rm = TRUE))
    expect_equivalent(coef2(eSSC.res), GS, tol = 1e-6)
    ## sigma(e.GS)^2
    expect_equivalent(vcov2(eSSC.lvm)[1:2,1:2], vcov(e.GS), tol = 1e-6)

    expect_equal(vcov2(eSSC.GS)[1:2,1:2],vcov(e.GS)[1:2,1:2], tol = 1e-6)
    expect_equal(vcov2(eSSC.lm1)[1,1],vcov(e.lm1)[1,1], tol = 1e-6)
    expect_equal(vcov2(eSSC.lm2)[1,1],vcov(e.lm2)[1,1], tol = 1e-6)
    expect_equal(vcov2(eSSC.lvm),vcov(e.GS),tol = 1e-6)
})



## * Brice, 01/27/20 6:12, ssc residuals under constrains
## path <- "C:/Users/hpl802/Downloads"
## butils.base:::sourcePackage("lavaSearch2", path = path, c.code = TRUE, trace = TRUE) ## version 1.5.5

dd <- data.frame("Y1" = c(-0.35, -0.87, -2.24, -0.7, 0.04, -1.46, -1.29, 0.6, -1.44, -1.64, -0.33, 1.12, -2, 0.66, 0.09, 1.18, -1.72, -1.02, 1.76, -0.48, -0.63, -1.95, -0.98, -2.8, -0.61), 
                 "eta" = c(-0.37, -0.69, -0.87, -0.1, -0.25, -1.85, -0.08, 0.97, 0.18, -1.38, -1.44, 0.36, -1.76, -0.32, -0.65, 1.09, -0.76, -0.83, 0.83, -0.97, -0.03, 0.23, -0.3, -0.68, 0.66), 
                 "Y2" = c(-0.77, -1.02, 0.5, 2.04, 0.25, -1.07, -0.98, 1.5, -0.46, -1.09, -2.67, -0.09, -2.59, 0.02, 0.41, 2.3, -0.03, -1.31, 1.4, -2.21, 0.35, -1.2, -1.35, -0.9, -0.83))

m <- lvm(Y1[0:sigma1] ~ 1*eta,
         Y2[0:sigma2] ~ 1*eta,
         eta[mu:1]  ~ 1
         )
latent(m) <- ~eta

e <- estimate(m, dd)

test_that("bug in version 1.5.4 (incorrect handling of the constrain when computing Omega)", {
    
    test.res1 <- estimate2(e, df = "Satterthwaite", ssc = "residuals", derivative = "analytic")
    ## test.cox1 <- estimate2(e, df = "Satterthwaite", ssc = "Cox", derivative = "analytic")

    expect_equal(test.res1$sCorrect$param, c("eta" = -0.58362569, "Y1~~Y1" = 0.54761303, "Y2~~Y2" = 1.01146135),
                 tol = 1e-6)
    expect_equal(test.cox1$sCorrect$param, c("eta" = -0.58362569, "Y1~~Y1" = 0.50602817, "Y2~~Y2" = 0.95769503),
                 tol = 1e-6)
    
})

m <- lvm(Y1[mu:sigma1] ~ 1,
         Y2[mu:sigma2] ~ 1)
e <- estimate(m, dd)

## test_that("bug in version 1.5.4 (incorrect handling of the constrain when computing Omega)", {
    
##     test.res1 <- sCorrect(e, df = "Satterthwaite", ssc = "residuals", derivative = "analytic")
##     test.res2 <- sCorrect(e, df = "Satterthwaite", ssc = "residuals", derivative = "numeric")
##     summary2(test.res1)$table2
##     summary2(test.res2)$table2
##     range( test.res1$sCorrect$dInformation - test.res1$sCorrect$dInformation )
    
    
## })

## * Brice 04/15/20 9:26 multcomp
mSim <- lvm(Y1~X1,Y2~X1)
n <- 25

set.seed(10)
d <- lava::sim(mSim, n)
dNA <- rbind(d,c(Y1=NA,X1=1,Y2=2))
dNA$id <- 1:26
ls.lmALL <- list("Y1" = estimate(lvm(Y1~X1), data = d),
                 "Y2" = estimate(lvm(Y2~X1), data = d))
ls.lmRED <- list("Y1" = estimate(lvm(Y1~X1), data = d))
ls.lmNA <- list("Y1" = estimate(lvm(Y1~X1), data = dNA),
                 "Y2" = estimate(lvm(Y2~X1), data = dNA))
class(ls.lmALL) <- "mmm"
class(ls.lmRED) <- "mmm"
class(ls.lmNA) <- "mmm"


test_that("Stability by subset", {
    glht.ALL <- glht2(ls.lmALL, linfct = "X1=0")
    glht.RED <- glht2(ls.lmRED, linfct = "X1=0")
    glht.NA <- glht2(ls.lmNA, linfct = "X1=0", cluster = "id")

    index.model1 <- which(grepl("Y1: ",colnames(glht.ALL$vcov)))
    expect_equal(glht.ALL$vcov[index.model1,index.model1], glht.RED$vcov, tol = 1e-6)
    expect_equal(glht.RED$vcov[index.model1,index.model1], glht.NA$vcov[index.model1,index.model1], tol = 1e-6)
})

######################################################################
### test-previousBug.R ends here

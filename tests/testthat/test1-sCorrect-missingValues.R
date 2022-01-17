### test1-sCorrect-missingValues.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: mar  7 2018 (13:39) 
## Version: 
## Last-Updated: jan 17 2022 (17:08) 
##           By: Brice Ozenne
##     Update #: 51
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
## rm(list = ls())
if(FALSE){ ## already called in test-all.R
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))
library(nlme)
context("sCorrect (dealing with missing values)")

## * simulation
n <- 2e1
mSim <- lvm(c(Y1~eta1,Y2~eta1+X2,Y3~eta1+X1,
           Z1~eta2,Z2~eta2,Z3~eta2+X3))
regression(mSim) <- eta1~X1+Gender
latent(mSim) <- ~eta1+eta2
categorical(mSim, labels = c("Male","Female")) <- ~Gender
transform(mSim, Id~Y1) <- function(x){1:NROW(x)}
set.seed(10)
d <- lava::sim(mSim, n = n, latent = FALSE)
dL <- reshape2::melt(d, id.vars = c("Id","X1","X2","X3","Gender"),
           measure.vars = c("Y1","Y2","Y3","Z1","Z2","Z3"))
dLred <- dL[dL$variable %in% c("Y1","Y2","Y3"),]
dLred2 <- dL[dL$variable %in% c("Y1","Y2"),]
dLred3 <- dLred2
dLred3[dL$variable == "Y2","Id"] <- paste0("-",dLred3[dL$variable == "Y2","Id"])

## ** t test
## formula:
## df = \frac{ 2 * s_pool^2 }{ var(s_pool^2) }
##    = \frac{ ( s_X^2/m + s_Y^2/n )^2}{( s_X^4/(m(m-1)) + s_Y^4/(n(n-1)))}

## using the t test function
e.ttest <- t.test(d$Y1, d$Y2)
e.ttest$parameter

## by hand
sX1 <- var(d$Y1)/n
sX2 <- var(d$Y2)/n
df <- (sX1+sX2)^2/(sX1^2/(n-1) + sX2^2/(n-1))
df-e.ttest$parameter

## * LVM: factor model
m <- lvm(c(Y1~eta1,Y2~eta1,Y3~eta1+X1))
regression(m) <- eta1~X1+X2

e.lvm <- estimate(m,d)
e2.lvm <- estimate2(e.lvm)

## ** complete case analysis
missing.Row <- d[1,]
missing.Row[,"Id"] <- -1
missing.Row[,c("Y1","Y2","Y3")] <- NA
## eNA.lvm <- estimate(m, rbind(d,missing.Row), missing = FALSE)
eNA.lvm <- estimate(m, rbind(missing.Row,d,missing.Row), missing = FALSE)

test_that("complete case analysis (factor model)", {
    eNA2.lvm <- estimate2(eNA.lvm)
    
    expect_equal(model.tables(eNA2.lvm), model.tables(e2.lvm))
    ## estimate        se        df       lower     upper   statistic     p.value
    ## eta1       -0.37387088 0.2962948 18.338440 -0.99554076 0.2477990 -1.26182045 0.222827159
    ## Y2         -0.02252887 0.3297457 15.258888 -0.72432797 0.6792702 -0.06832194 0.946416628
    ## Y3          0.38272845 0.2851559  5.682026 -0.32459821 1.0900551  1.34217267 0.230678148
    ## eta1~X1     0.99599616 0.3134807 18.469394  0.33859531 1.6533970  3.17721651 0.005092730
    ## eta1~X2    -0.04275890 0.2607688  9.859874 -0.62490927 0.5393915 -0.16397243 0.873065310
    ## Y2~eta1     1.05707590 0.2723211  5.481730  0.37515670 1.7389951  3.88172564 0.009722483
    ## Y3~eta1     1.08664682 0.3566308  1.982418 -0.46092704 2.6342207  3.04697992 0.093951240
    ## Y3~X1       0.61495545 0.4296209  3.045199 -0.74088088 1.9707918  1.43139085 0.246430988
    ## Y1~~Y1      0.49861889 0.3398654  5.101312 -0.36983983 1.3670776          NA          NA
    ## eta1~~eta1  1.10299240 0.5611968  4.190452 -0.42760479 2.6335896          NA          NA
    ## Y2~~Y2      1.60214650 0.6249695  4.823447 -0.02223785 3.2265309          NA          NA
    ## Y3~~Y3      0.64437632 0.4273389  4.023035 -0.53943230 1.8281849          NA          NA
})

## ** full information
missing.Row <- d[1,]
missing.Row[,"Id"] <- -1
missing.Row[,c("Y1","Y2")] <- NA
eNA.lvm <- estimate(m, rbind(missing.Row,d,missing.Row), missing = TRUE)


test_that("full information (factor model)", {
    expect_equal(information(eNA.lvm), unname(information2(eNA.lvm, ssc = FALSE)), tol = 1e-6)

    ## NOT TRUE!!!! issue in lava or on purpose?
    ## expect_equal(unname(information(eNA.lvm)), unname(solve(vcov(eNA.lvm))), tol = 1e-6)

    eNA2.lvm <- estimate2(eNA.lvm)
    model.tables(eNA2.lvm)
    ##               estimate        se        df      lower     upper  statistic     p.value
    ## eta1       -0.38547381 0.2851540 19.464112 -0.9813462 0.2103985 -1.3518095 0.191937640
    ## Y3          0.35963220 0.2808357  4.987241 -0.3628347 1.0820991  1.2805789 0.256661658
    ## Y2         -0.02211643 0.3302074 15.101848 -0.7255237 0.6812908 -0.0669774 0.947478382
    ## eta1~X1     1.01140559 0.2987973 20.052756  0.3882305 1.6345807  3.3849225 0.002933646
    ## eta1~X2    -0.02504732 0.2438504  9.187004 -0.5749683 0.5248737 -0.1027159 0.920395648
    ## Y3~eta1     1.05020734 0.3633959  1.961287 -0.5433094 2.6437241  2.8899816 0.104091914
    ## Y3~X1       0.66179756 0.4338628  2.801850 -0.7758670 2.0994621  1.5253613 0.230806721
    ## eta1~~eta1  1.02653446 0.5313942  3.405538 -0.5561885 2.6092574         NA          NA
    ## Y3~~Y3      0.61722688 0.4035118  3.918306 -0.5123737 1.7468274         NA          NA
    ## Y2~eta1     1.05999331 0.2760511  5.764207  0.3777659 1.7422208  3.8398445 0.009242974
    ## Y1~~Y1      0.50327981 0.3397094  4.999472 -0.3699987 1.3765584         NA          NA
    ## Y2~~Y2      1.59766033 0.6235258  4.836780 -0.0215698 3.2168905         NA          NA
    expect_equal(summary(eNA2.lvm)$table2$df,
                 c(19.46411161, 4.98724128, 15.10184787, 20.05275552, 9.18700366, 1.96128742, 2.80184994, 3.4055378, 3.91830642, 5.7642073, 4.99947214, 4.83678),
                 tol = 1e-6)
})

##----------------------------------------------------------------------
### test1-sCorrect-missingValues.R ends here

### test-previousBug.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 19 2019 (10:17) 
## Version: 
## Last-Updated: jan 28 2020 (18:06) 
##           By: Brice Ozenne
##     Update #: 18
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * Brice, 11/19/19 10:18, summary2 (name coef)
long <- data.frame("id" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"), 
                   "visit" = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4), 
                   "weight" = c(127.2, 165.2, 109.7, 146.2, 113.1, 158.8, 115.4, 123.8, 105.8, 118.0, 137.7, 129.0, 117.4, 118.0, 115.0, 131.2, 122.4, 151.6, 173.0, 100.9, 120.7, 153.4, 101.6, 142.4, 105.6, 143.6, 108.8, 115.1, 102.0, 105.3, 131.9, 120.7, 108.5, 113.4, 109.7, 127.3, 113.9, 143.0, 162.2,  95.7, 115.5, 149.2,  97.7, 136.7,  99.9, 134.6, 103.0, 112.2,  98.0,  99.3, 126.3, 116.6, 104.5, 108.3, 103.5, 119.5, 109.0, 135.3, 155.0,  89.9, 108.1, 132.0,  87.1, 123.0,  87.7, 108.7,  88.9,  97.6,  86.6,  90.9, 105.4, 100.0,  95.4,  93.6,  94.1, 103.5,  99.4, 118.5, 148.0,  78.8), 
                   "glucagonAUC" = c( 5032.50, 12142.50, 10321.35,  6693.00,  7090.50, 10386.00, 10258.80, 16797.75,  2500.50, 11862.00,  5199.00,  5146.50,  3642.00,  4372.95,  5410.50,  7015.50,  6673.50,  5485.50,  6879.00, 14299.50,  4942.50, 14083.50,  6202.50,  6631.50, NA,  7609.50,  8424.60, 16300.50,  2376.00,  8212.50,  5502.00,  5217.00,  5911.50,  4520.10,  7833.00,  5497.50,  4857.00,  5010.00,  7953.00,  8739.00, 20421.00, 10945.50, 20121.00, 13090.50, 19155.00, 11778.00, 29979.75, 11040.00, 13657.50, 22875.00,  7906.50, 12897.00, 26555.10, 12903.90, NA, 16290.00, 17560.50, 16269.00, 12036.00, 26638.50,  9249.45,  7612.50, 17704.50,  4551.00, 12345.00,  8014.80, 11837.70,  6163.50, 11286.00, 12339.00,  6543.00, 11499.00, 23245.50, 12536.10, 18148.50, 10536.00,  8434.95,  7441.50, 10362.00, 11410.50), 
                   "time" = c("-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month"), 
                   "time2" = c("-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-3 month", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "-1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+1 week", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month", "+3 month"))

## long$time3 <- as.factor(as.numeric(long$time))
fit.time <- gls(weight~time3,
                data=long,
                correlation=corSymm(form=~visit|id),
                weights=varIdent(form=~1|time),
                na.action=na.exclude,
                control=glsControl(opt="optim"))
summary2(fit.time)


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

ddL <- dcast(ddW, value.var = "neocortex.log", id~group)
names(ddL) <- c("id","G1","G2")
m <- lvm(G1~1,G2~1)
e.lvm <- estimate(m, ddL)

test_that("sCorrect in stratified GLS equivalent to separate LM", {
    eSSC.GS <- sCorrect(e.GS, ssc = "Cox")
    eSSC0.GS <- sCorrect(e.GS, ssc = "residuals")
    
    eSSC.lm1 <- sCorrect(e.lm1, ssc = "Cox")
    eSSC.lm2 <- sCorrect(e.lm2, ssc = "Cox")
    eSSC.lvm <- sCorrect(e.lvm, ssc = "Cox")

    GS <- c(mu1 = mean(ddL$G1, na.rm = TRUE),
            mu2 = mean(ddL$G2, na.rm = TRUE),
            sigma1 = var(ddL$G1, na.rm = TRUE),
            sigma2 = var(ddL$G2, na.rm = TRUE))
    expect_equivalent(coef2(eSSC.lvm), GS, tol = 1e-6)
    ## sigma(e.GS)^2
    expect_equivalent(vcov2(eSSC.lvm)[1:2,1:2], vcov(e.GS), tol = 1e-6)

    eSSC.GS$sCorrect$ssc$param0
    diag(eSSC.GS$sCorrect$ssc$Psi)

    eSSC.lm1$sCorrect$ssc$param0
    eSSC.lm1$sCorrect$ssc$Psi
    eSSC.lm1$sCorrect$param

    eSSC.lvm$sCorrect$ssc$param0
    coef2(eSSC.GS)
    coef2(eSSC.lvm)

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

test_that("Comparing the two corrections (residuals,Cox): simular but not equal values", {
    test <- sCorrect(e, df = "Satterthwaite", ssc = "residuals", derivative = "analytic")
    test2 <- sCorrect(e, df = "Satterthwaite", ssc = "Cox", derivative = "analytic")

    expect_equal(test$sCorrect$param, c("eta" = -0.58362569, "Y1~~Y1" = 0.54761303, "Y2~~Y2" = 1.01146135),
                 tol = 1e-6)
    expect_equal(test2$sCorrect$param, c("eta" = -0.58362569, "Y1~~Y1" = 0.50602817, "Y2~~Y2" = 0.95769503),
                 tol = 1e-6)
    
})


######################################################################
### test-previousBug.R ends here

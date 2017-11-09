### test-coefType.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (14:52) 
## Version: 
## last-updated: nov  9 2017 (16:32) 
##           By: Brice Ozenne
##     Update #: 56
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

library(testthat)

context("coefType")

lava.options(symbols = c("~","~~"))

## * linear regression
## ** only continous variables
m <- lvm(Y~X1+X2)
e <- estimate(m, sim(m, 1e2))


dt.truth <- data.table(name = c("Y","Y~X1","Y~X2","Y~~Y"),
                       type = c("intercept","regression","regression","variance"),
                       letter = c("nu","K","K","Sigma_var")
                       )

test_that("coefType - lm", {

    test <- coefType(m, as.lava = TRUE)
    GS <- setNames(dt.truth$type,dt.truth$name)
    expect_equal(test,GS)

    test <- coefType(m, as.lava = FALSE)
    expect_equal(test[!is.na(lava),detail],dt.truth$letter)
    expect_equal(test[!is.na(lava),name],dt.truth$name)

    test <- coefType(e, as.lava = TRUE)
    GS <- setNames(dt.truth$type,dt.truth$name)[names(test)]
    expect_equal(test,GS)

    test <- coefType(e, as.lava = FALSE)
    expect_equal(test[!is.na(lava),detail],dt.truth$letter)
    expect_equal(test[!is.na(lava),name],dt.truth$name)

})

test_that("coefCov - lm", {
    expect_true(is.null(coefCov(m)))
    expect_true(is.null(coefCov(m, value = TRUE)))

    expect_equal(coefCov(m, keep.var = TRUE),4)
    expect_equal(coefCov(m, value = TRUE, keep.var = TRUE),"Y~~Y")
})

test_that("coefIndexModel - lm", {
    expect_true(all(coefIndexModel(m)==1))
    expect_true(all(coefIndexModel(e)==1))
})

test_that("coefIntercept - lm", {
    expect_equal(coefIntercept(m, keep.var = TRUE),1)
    expect_equal(coefIntercept(m, value = TRUE, keep.var = TRUE),"Y")
}) 

test_that("coefReg - lm", {
    expect_equal(coefReg(m, keep.var = TRUE),2:3)
    expect_equal(coefReg(m, value = TRUE, keep.var = TRUE),c("Y~X1","Y~X2"))
}) 

## ** 1 categorical variable
m <- lvm(Y~X1+X2+X3)
mSim <- m
categorical(mSim, labels =c("a","b","c")) <- ~X2
d <- sim(mSim, 1e2)
e <- estimate(m, d)

dt.truth <- data.table(name = c("Y","Y~X1","Y~X2b","Y~X2c","Y~X3","Y~~Y"),
                       type = c("intercept","regression","regression","regression","regression","variance"),
                       letter = c("nu","K","K","K","K","Sigma_var")
                       )

test_that("coefType - lm", {
    
    test <- coefType(mSim, as.lava = TRUE)
    expect_equal(names(test),as.character(coef(mSim)))

    test <- coefType(e, as.lava = TRUE)
    GS <- setNames(dt.truth$type,dt.truth$name)[names(coef(e))]
    expect_equal(test,GS)
    
    test <- coefType(e, as.lava = FALSE)
    expect_equal(test[!is.na(lava),detail],dt.truth$letter)
    expect_equal(test[!is.na(lava),name],dt.truth$name)
   
})

## ** Several categorical variables
m <- lvm(Y1~X1+X2,Y2~X1+X3)
mSim <- m
categorical(mSim, labels = c("a","b","c")) <- "X1"
d <- sim(mSim, 1e2)
e <- estimate(m, d)

dt.truth <- data.table(name = c("Y1","Y2",
                                "Y1~X1b","Y1~X1c","Y1~X2",
                                "Y2~X1b","Y2~X1c","Y2~X3",
                                "Y1~~Y1","Y2~~Y2"))
dt.truth[,type := c(rep("intercept",2),rep("regression",6),rep("variance",2))]
dt.truth[,letter := c(rep("nu",2),rep("K",6),rep("Sigma_var",2))]

test_that("coefType - lm", {

    test <- coefType(mSim, as.lava = TRUE)
    expect_equal(names(test),as.character(coef(mSim)))

    test <- coefType(e, as.lava = TRUE)
    GS <- setNames(dt.truth$type,dt.truth$name)[names(coef(e))]
    expect_equal(test,GS)
    
    test <- coefType(e, as.lava = FALSE)
    expect_equal(test[!is.na(lava),detail],dt.truth$letter)
    expect_equal(test[!is.na(lava),name],dt.truth$name)
})

## * lvm
## ** categorical variable
m <- lvm(c(Y1,Y2,Y3)~eta1, c(Z1,Z2,Z3)~eta2)
regression(m) <- Y1~X1
regression(m) <- eta2~eta1+X2
covariance(m) <- Y1~Z1
latent(m) <- ~eta1+eta2
mSim <- m
categorical(mSim, labels =c("a","b","c")) <- ~X2
d <- sim(mSim, 1e2, latent = FALSE)
e <- estimate(m, d)

dt.truth <- rbind(c("Y2","intercept","nu"),
                  c("Y3","intercept","nu"),
                  c("Z2","intercept","nu"),
                  c("Z3","intercept","nu"),
                  ##
                  c("eta1","intercept","alpha"),
                  c("eta2","intercept","alpha"),
                  ##
                  c("Y1~X1","regression","K"),
                  ##
                  c("Y2~eta1","regression","Lambda"),
                  c("Y3~eta1","regression","Lambda"),
                  c("Z2~eta2","regression","Lambda"),
                  c("Z3~eta2","regression","Lambda"),
                  ##
                  c("eta2~eta1","regression","B"),
                  ##
                  c("eta2~X2b","regression","Gamma"),
                  c("eta2~X2c","regression","Gamma"),
                  ##
                  c("Y1~~Y1","variance","Sigma_var"),
                  c("Y2~~Y2","variance","Sigma_var"),
                  c("Y3~~Y3","variance","Sigma_var"),
                  c("Z1~~Z1","variance","Sigma_var"),
                  c("Z2~~Z2","variance","Sigma_var"),
                  c("Z3~~Z3","variance","Sigma_var"),
                  ##
                  c("eta1~~eta1","variance","Psi_var"),
                  c("eta2~~eta2","variance","Psi_var"),
                  ##
                  c("Y1~~Z1","covariance","Sigma_cov")
                  )


dt.truth <- as.data.table(dt.truth)
names(dt.truth) <- c("name","type","letter")
setkeyv(dt.truth, c("type","letter","name"))

test_that("coefType - lvm", {

    test <- coefType(mSim, as.lava = TRUE)
    expect_equal(names(test),as.character(coef(mSim)))
    
    test <- coefType(e, as.lava = TRUE)
    GS <- setNames(dt.truth$type,dt.truth$name)[names(coef(e))]
    expect_equal(test,GS)

    test <- coefType(e, as.lava = FALSE)
    expect_equal(test[!is.na(lava),detail],dt.truth$letter)
    expect_equal(test[!is.na(lava),name],dt.truth$name)
})

## ** constrains (0 mean 1 loading)

m <- lvm(c(Y1~0+1*eta1,Y2~0+1*eta1,Y3~0+1*eta1,
           Z1~0+1*eta2,Z2~0+1*eta2,Z3~0+1*eta2))
latent(m) <- ~eta1 + eta2
covariance(m) <- Y1~Y2

e <- estimate(m, sim(m,1e2))

dt.truth <- rbind(data.table("Y1","intercept","nu",TRUE),
                  data.table("Y2","intercept","nu",TRUE),
                  data.table("Y3","intercept","nu",TRUE),
                  data.table("Z1","intercept","nu",TRUE),
                  data.table("Z2","intercept","nu",TRUE),
                  data.table("Z3","intercept","nu",TRUE),
                  ##
                  data.table("eta1","intercept","alpha",FALSE),
                  data.table("eta2","intercept","alpha",FALSE),
                  ##
                  data.table("Y1~eta1","regression","Lambda",TRUE),
                  data.table("Y2~eta1","regression","Lambda",TRUE),
                  data.table("Y3~eta1","regression","Lambda",TRUE),
                  data.table("Z1~eta2","regression","Lambda",TRUE),
                  data.table("Z2~eta2","regression","Lambda",TRUE),
                  data.table("Z3~eta2","regression","Lambda",TRUE),
                  ##
                  data.table("Y1~~Y1","variance","Sigma_var",FALSE),
                  data.table("Y2~~Y2","variance","Sigma_var",FALSE),
                  data.table("Y3~~Y3","variance","Sigma_var",FALSE),
                  data.table("Z1~~Z1","variance","Sigma_var",FALSE),
                  data.table("Z2~~Z2","variance","Sigma_var",FALSE),
                  data.table("Z3~~Z3","variance","Sigma_var",FALSE),
                  ##
                  data.table("eta1~~eta1","variance","Psi_var",FALSE),
                  data.table("eta2~~eta2","variance","Psi_var",FALSE),
                  ##
                  data.table("Y1~~Y2","covariance","Sigma_cov",FALSE)
                  )
dt.truth <- as.data.table(dt.truth)
names(dt.truth) <- c("name","type","letter","fixed")
setkeyv(dt.truth, c("type","letter","name"))

test_that("coefType - constrains 0/1", {

    test <- coefType(m, as.lava = TRUE)
    expect_equal(names(test),as.character(coef(m)))

    test <- coefType(e, as.lava = TRUE)
    GS <- setNames(dt.truth[fixed==FALSE,type],dt.truth[fixed==FALSE,name])[names(coef(e))]
    expect_equal(test,GS)
    
    test <- coefType(e, as.lava = FALSE)
    expect_equal(test$detail,dt.truth$letter)
    expect_equal(test$name,dt.truth$name)    
})


## ** constrains (common mean and variance)

m <- lvm(c(Y1[mu:sigma]~beta*eta1,Y2[mu:sigma]~1*eta1,Y3~1*eta1+X1),eta1~X1)
mSim <- m
categorical(mSim, labels = as.character(1:3)) <- ~X1

d <- sim(mSim, latent=FALSE, 1e2)
e <- estimate(m, d)

test_that("coefType - constrains mean/variance", {
    test <- coefType(m, as.lava = TRUE)
    expect_equal(names(test),as.character(coef(m)))

    test <- coefType(e, as.lava = TRUE)
    #GS <- setNames(dt.truth[fixed==FALSE,type],dt.truth[fixed==FALSE,name])[names(coef(e))]
    #expect_equal(test,GS)
 
    coefType(e, as.lava = FALSE)
})

 
#----------------------------------------------------------------------
### test-coefType.R ends here


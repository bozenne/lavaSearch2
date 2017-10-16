### test-coefType.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 12 2017 (14:52) 
## Version: 
## last-updated: okt 16 2017 (16:31) 
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

context("coefType")

## * linear regression

## ** only continous variables
m <- lvm(Y~X1+X2)
e <- estimate(m, sim(m, 1e2))


dt.truth <- data.table(name = c("Y","Y~X1","Y~X2","Y~~Y"),
                       type = c("intercept","regression","regression","variance"),
                       letter = c("nu","K","K","Sigma_var")
                       )

test_that("coefType - lm", {

    test <- coefType(m, detailed = FALSE)
    GS <- setNames(dt.truth$type,dt.truth$name)
    expect_equal(test,GS)

    test <- coefType(m, detailed = TRUE)
    GS <- setNames(dt.truth$letter,dt.truth$name)
    expect_equal(test,GS)

    test <- coefType(e, detailed = FALSE)
    GS <- setNames(dt.truth$type,dt.truth$name)[names(test)]
    attr(GS,"reference") <- setNames(rep(FALSE,length(GS)),names(GS))
    expect_equal(test,GS)

    test <- coefType(e, detailed = TRUE)
    GS <- setNames(dt.truth$letter,dt.truth$name)[names(test)]
    attr(GS,"reference") <- setNames(rep(FALSE,length(GS)),names(GS))
    expect_equal(test,GS)
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

## ** categorical variable
m <- lvm(Y~X1+X2+X3)
m.Sim <- m
categorical(m.Sim, labels =c("a","b","c")) <- ~X2
d <- sim(m.Sim, 1e2)
e <- estimate(m, d)


dt.truth <- data.table(name = c("Y","Y~X1","Y~X3","Y~X2b","Y~X2c","Y~~Y"),
                       type = c("intercept","regression","regression","regression","regression","variance"),
                       letter = c("nu","K","K","K","K","Sigma_var")
                       )

test_that("coefType - lm", {

    test <- coefType(e, detailed = FALSE)
    GS <- setNames(dt.truth$type,dt.truth$name)
    expect_equal(length(test),length(GS))
    GS <- GS[names(test)]
    attr(GS,"reference") <- setNames(rep(FALSE,length(GS)),names(GS))
    expect_equal(test,GS)
    
    test <- coefType(e, detailed = TRUE)
    GS <- setNames(dt.truth$letter,dt.truth$name)
    expect_equal(length(test),length(GS))
    GS <- GS[names(test)]
    attr(GS,"reference") <- setNames(rep(FALSE,length(GS)),names(GS))
    expect_equal(test,GS)
    
})

## * lvm

## ** categorical variable
m <- lvm(c(Y1,Y2,Y3)~eta1, c(Z1,Z2,Z3)~eta2)
regression(m) <- Y1~X1
regression(m) <- eta2~eta1+X2
covariance(m) <- Y1~Z1
latent(m) <- ~eta1+eta2
m.Sim <- m
categorical(m.Sim, labels =c("a","b","c")) <- ~X2
d <- sim(m.Sim, 1e2, latent = FALSE)
e <- estimate(m, d)

dt.truth <- rbind(c("Y1","intercept","nu"),
                  c("Y2","intercept","nu"),
                  c("Y3","intercept","nu"),
                  c("Z1","intercept","nu"),
                  c("Z2","intercept","nu"),
                  c("Z3","intercept","nu"),
                  c("eta1","intercept","alpha"),
                  c("eta2","intercept","alpha"),
                  ##
                  c("Y1~eta1","regression","Lambda"),
                  c("Y1~X1","regression","K"),
                  c("Y2~eta1","regression","Lambda"),
                  c("Y3~eta1","regression","Lambda"),
                  c("Z1~eta2","regression","Lambda"),
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
                  c("Y1~~Z1","covariance","Sigma_cov"),
                  ##
                  c("eta2~~eta2","variance","Psi_var"),
                  c("eta1~~eta1","variance","Psi_var")
                  )


dt.truth <- as.data.table(dt.truth)
names(dt.truth) <- c("name","type","letter")

test_that("coefType - lm", {

    test <- coefType(e, detailed = FALSE)
    GS <- setNames(dt.truth$type,dt.truth$name)
    expect_equal(length(test),length(GS))
    GS <- GS[names(test)]
    attr(GS,"reference") <- setNames(rep(FALSE,length(GS)),names(GS))
    attr(GS,"reference")["Y1"] <- TRUE
    attr(GS,"reference")["Y1~eta1"] <- TRUE
    attr(GS,"reference")["Z1"] <- TRUE
    attr(GS,"reference")["Z1~eta2"] <- TRUE
    expect_equal(test,GS)
    
    test <- coefType(e, detailed = TRUE)
    GS <- setNames(dt.truth$letter,dt.truth$name)
    expect_equal(length(test),length(GS))
    GS <- GS[names(test)]
    attr(GS,"reference") <- setNames(rep(FALSE,length(GS)),names(GS))
    attr(GS,"reference")["Y1"] <- TRUE
    attr(GS,"reference")["Y1~eta1"] <- TRUE
    attr(GS,"reference")["Z1"] <- TRUE
    attr(GS,"reference")["Z1~eta2"] <- TRUE
    expect_equal(test,GS)
    
})

#----------------------------------------------------------------------
### test-coefType.R ends here

## ** constrains

m <- lvm(c(Y1~0+1*eta1,Y2~0+1*eta1,Y3~0+1*eta1,
           Z1~0+1*eta2,Z2~0+1*eta2,Z3~0+1*eta2))
latent(m) <- ~eta1 + eta2
#covariance(m,var1="Y1",var2="Y2") <- 0.5
covariance(m) <- Y1~Y2

e <- estimate(m, sim(m,1e2))
coefType(m)
coef(e, level = 9)
coefType(e)

dt.truth <- rbind(c("Y1","intercept","nu"),
                  c("Y2","intercept","nu"),
                  c("Y3","intercept","nu"),
                  c("Z1","intercept","nu"),
                  c("Z2","intercept","nu"),
                  c("Z3","intercept","nu"),
                  c("eta1","intercept","alpha"),
                  c("eta2","intercept","alpha"),
                  ##
                  c("Y1~eta1","regression","Lambda"),
                  c("Y2~eta1","regression","Lambda"),
                  c("Y3~eta1","regression","Lambda"),
                  c("Z1~eta2","regression","Lambda"),
                  c("Z2~eta2","regression","Lambda"),
                  c("Z3~eta2","regression","Lambda"),
                  ##
                  c("Y1~~Y1","variance","Sigma_var"),
                  c("Y2~~Y2","variance","Sigma_var"),
                  c("Y3~~Y3","variance","Sigma_var"),
                  c("Z1~~Z1","variance","Sigma_var"),
                  c("Z2~~Z2","variance","Sigma_var"),
                  c("Z3~~Z3","variance","Sigma_var"),
                  ##
                  c("Y1~~Y2","covariance","Sigma_cov"),
                  ##
                  c("eta2~~eta2","variance","Psi_var"),
                  c("eta1~~eta1","variance","Psi_var")
                  )
dt.truth <- as.data.table(dt.truth)
names(dt.truth) <- c("name","type","letter")


test_that("coefType - constrains", {

    test <- coefType(e, detailed = FALSE)
    GS <- setNames(dt.truth$type,dt.truth$name)
    expect_equal(length(test),length(GS))
    GS <- GS[names(test)]
    attr(test,"reference") <- NULL
    expect_equal(test,GS)
    
    test <- coefType(e, detailed = TRUE)
    GS <- setNames(dt.truth$letter,dt.truth$name)
    expect_equal(length(test),length(GS))
    GS <- GS[names(test)]
    attr(test,"reference") <- NULL
    expect_equal(test,GS)
    
})

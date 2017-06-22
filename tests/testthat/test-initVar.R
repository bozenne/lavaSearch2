library(testthat)
# library(lavaReduce)
# library(lava)

context("#### initVar #### \n")

initVar_link <-  lavaReduce:::initVar_link

lava.options(symbols = c("~","~~"))

test_that("initVar - two different variables",{
  l <- list(var1 = "a", var2 = "b")  
  f <- a~b
  # regression
  expect_equal(l, initVar_link(var1 = a~b))
  expect_equal(l, initVar_link(var1 = a ~ b))
  expect_equal(list(f),initVar_link(var1 = "a ~ b", format = "formula"))
  # covariance
  expect_equal(l, initVar_link(var1 = "a~~b"))
  # var
  expect_equal(l, initVar_link(var1 = "a", var2 = "b"))
})

test_that("initVar - one repeated variable",{
  
  l <- list(var1 = "X1", var2 = "X1")
  # regression
  expect_equal(l, initVar_link(var1 = X1~X1))
  expect_equal(l, initVar_link(var1 = "X1~~X1"))
})

test_that("initVar - no response variable",{
  
  l1 <- list(var1 = NULL, var2 = "X1")
  l12 <- list(var1 = NULL, var2 = c("X1","X2"))
  # regression
  expect_equal(l1, initVar_link(var1 = ~X1))
  expect_equal(l12, initVar_link(var1 = ~X1+X2))
})


initVar_link(var1 = a ~ b+c+d*e, format = "list")
initVar_link(var1 = a ~ b+c+d*e, format = "txt.formula")
initVar_link(var1 = a ~ b+c+d*e, format = "formula")




initVar_link(var1 = Y~X1+X2)
initVar_link(var1 = Y~X1+X2, repVar1 = TRUE)
initVar_link(var1 = Y~X1+X2, format = "formula")
initVar_link(var1 = Y~X1+X2, format = "txt.formula")

lava.options(symbols = c("<-","<->"))
initVar_link(var1 = "Y<-X1+X2", repVar1 = TRUE)
initVar_link(var1 = "Y<-X1+X2", format = "formula")
initVar_link(var1 = "Y<-X1+X2", format = "txt.formula")




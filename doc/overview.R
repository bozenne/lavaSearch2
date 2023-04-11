## chunk 2
library(lavaSearch2)

## * Inference
## ** Introductory example

## chunk 3
## simulate data
mSim <- lvm(Y[1:1]~0.3*X1+0.2*X2)
set.seed(10)
df.data <- sim(mSim, 2e1)

## fit linear model
summary(lm(Y~X1+X2, data = df.data))$coef

## chunk 4
## fit latent variable model
m <- lvm(Y~X1+X2)
e <- estimate(m, data = df.data)

## extract Wald tests
summary(e)$coef

## chunk 5
summary2(e)$coef

## ** How it works in a nutshell
## ** Single univariate Wald test

## chunk 6
## simulate data
mSim <- lvm(Y1~eta,Y2~eta,Y3~0.4+0.4*eta,Y4~0.6+0.6*eta,eta~0.5*X1+0.7*X2)
latent(mSim) <- ~eta
set.seed(12)
df.data <- sim(mSim, n = 3e1, latent = FALSE)

## display
head(df.data)

## chunk 7
m <- lvm(c(Y1,Y2,Y3,Y4)~eta, eta~X1+X2)
e <- estimate(m, data = df.data)

## chunk 8
summary(e)$coef[c("Y2","Y3","Y2~eta","Y3~eta","eta~X1","eta~X2"), ]

## chunk 9
summary2(e)$coef[c("Y2","Y3","Y2~eta","Y3~eta","eta~X1","eta~X2"), ]

## chunk 10
summary2(e, ssc = FALSE)$coef[c("Y2","Y3","Y2~eta","Y3~eta","eta~X1","eta~X2"), ]

## chunk 11
summary2(e, df = FALSE)$coef[c("Y2","Y3","Y2~eta","Y3~eta","eta~X1","eta~X2"), ]

## ** Saving computation time with =estimate2=

## chunk 12
system.time(
    res <- summary2(e, ssc = FALSE)
)

## chunk 13
e2 <- estimate2(e)

## chunk 14
class(e2)

## chunk 15
system.time(
    summary(e2)
)

## ** Single multivariate Wald test

## chunk 16
resTest0 <- lava::compare(e, par = c("Y2","Y2~eta","eta~X1"))
resTest0

## chunk 17
resTest1 <- compare2(e, linfct = c("Y2","Y2~eta","eta~X1"))
resTest1

## chunk 18
resC <- createContrast(e, linfct = c("Y2=0","Y2~eta=0","eta~X1=0"))
resC$contrast

## chunk 19
resTest2 <- compare2(e2, linfct = resC$contrast)
identical(resTest1,resTest2)

## chunk 20
resTest3 <- compare2(e, ssc = FALSE, linfct = resC$contrast)
resTest3

## chunk 21
resTest0$statistic/qr(resC$contrast)$rank

## ** Robust Wald tests

## chunk 22
summary2(e, robust = TRUE)$coef[c("Y2","Y3","Y2~eta","Y3~eta","eta~X1","eta~X2"), ]

## chunk 23
compare2(e2, linfct = c("Y2","Y2~eta","eta~X1"), robust = TRUE)

## chunk 24
rbind(robust = diag(crossprod(iid(e))),
      model = diag(vcov(e)))

## ** Assessing the type 1 error of the testing procedure

## chunk 25
set.seed(10)
m.generative <- lvm(Y ~ X1 + X2 + Gene)
categorical(m.generative, labels = c("ss","ll")) <- ~Gene
d <- lava::sim(m.generative, n = 50, latent = FALSE)

## chunk 26
head(d)

## chunk 27
myModel <- lvm(Y ~ X1 + X2 + Gene)

## chunk 28
e <- estimate(myModel, data = d)
e

## chunk 29
mySimulation <- calibrateType1(e, 
                               param = "Y~Genell",
                               n.rep = 50, 
                               trace = FALSE, seed = 10)

## chunk 30
summary(mySimulation)

## * Adjustment for multiple comparisons
## ** Univariate Wald test, single model

## chunk 31
## simulate data
mSim <- lvm(Y ~ 0.25 * X1 + 0.3 * X2 + 0.35 * X3 + 0.4 * X4 + 0.45 * X5 + 0.5 * X6)
set.seed(10)
df.data <- sim(mSim, n = 4e1)

## fit lvm
e.lvm <- estimate(lvm(Y ~ X1 + X2 + X3 + X4 + X5 + X6), data = df.data)
name.coef <- names(coef(e.lvm))
n.coef <- length(name.coef)

## Create contrast matrix
resC <- createContrast(e.lvm, linfct = paste0("Y~X",1:6), rowname.rhs = FALSE)
resC$contrast

## chunk 32
e.glht <- multcomp::glht(e.lvm, linfct = resC$contrast, rhs = resC$null)
summary(e.glht)

## chunk 33
e.glht2 <- glht2(e.lvm, linfct = resC$contrast, rhs = resC$null)
summary(e.glht2)

## chunk 34
summary(e.glht2, test = multcomp::adjusted("free"))

## ** Univariate Wald test, multiple models

## chunk 35
mSim <- lvm(X ~ Age + 0.5*Treatment,
            Y ~ Gender + 0.25*Treatment,
            c(Z1,Z2,Z3) ~ eta, eta ~ 0.75*treatment,
            Age[40:5]~1)
latent(mSim) <- ~eta
categorical(mSim, labels = c("placebo","SSRI")) <- ~Treatment
categorical(mSim, labels = c("male","female")) <- ~Gender

n <- 5e1
set.seed(10)
df.data <- sim(mSim, n = n, latent = FALSE)
head(df.data)

## chunk 36
lvmX <- estimate(lvm(X ~ Age + Treatment), data = df.data)
lvmY <- estimate(lvm(Y ~ Gender + Treatment), data = df.data)
lvmZ <- estimate(lvm(c(Z1,Z2,Z3) ~ 1*eta, eta ~ -1 + Treatment), 
                 data = df.data)

## chunk 37
mmm.lvm <- multcomp::mmm(X = lvmX, Y = lvmY, Z = lvmZ)

## chunk 38
lvm.glht2 <- glht2(mmm.lvm, linfct = "TreatmentSSRI")
summary(lvm.glht2)

## chunk 39
summary(lvm.glht2, test = multcomp::adjusted("none"))

## * Model diagnostic
## ** Detection of local dependencies

## chunk 40
## simulate data
mSim <- lvm(c(y1,y2,y3)~u, u~x1+x2)
latent(mSim) <- ~u
covariance(mSim) <- y2~y3
transform(mSim, Id~u) <- function(x){1:NROW(x)}
set.seed(10)
df.data <- lava::sim(mSim, n = 125, latent = FALSE)
head(df.data)

## chunk 41
## fit model
m <- lvm(c(y1,y2,y3)~u, u~x1)
latent(m) <- ~u
addvar(m) <- ~x2 
e.lvm <- estimate(m, data = df.data)

## chunk 42
resScore <- modelsearch2(e.lvm, alpha = 0.1, trace = FALSE)
displayScore <- summary(resScore)

## chunk 43
resScore0 <- modelsearch(e.lvm, silent = TRUE)
c(statistic = sqrt(max(resScore0$test[,"Test Statistic"])), 
  p.value = min(resScore0$test[,"P-value"]))

## chunk 44
data.frame(link = displayScore$table[,"link"],
           none = displayScore$table[,"p.value"],
           bonferroni = displayScore$table[,"p.value"]*displayScore$table[1,"nTests"],
           max = displayScore$table[,"adjusted.p.value"])

## chunk 45
autoplot(resScore, param = "u~x1")

## chunk 46
resRed <- modelsearch2(e.lvm, link = c("y1~~y2","y1~~y3","y2~~y3"), trace = FALSE)
print(resRed)

## chunk 47
findNewLink(e.lvm$model, type = "covariance")$link

## ** Checking that the names of the variables in the model match those of the data

## chunk 48
## simulate data
set.seed(10)
df.data <- sim(lvm(Y~X1+X2), 1e2)

## fit model
mWrong <- lvm(Y ~ X + X2)
eWrong <- estimate(mWrong, data = df.data)

## chunk 49
checkData(mWrong, df.data)

## chunk 50
## simulate data
set.seed(10)
mSim <- lvm(c(Y1,Y2,Y3)~eta)
latent(mSim) <- ~eta
df.data <- sim(mSim, n = 1e2, latent = FALSE)

## fit model
m <- lvm(c(Y1,Y2,Y3)~eta)
checkData(m, data = df.data)

## chunk 51
latent(m) <- ~eta
checkData(m, data = df.data)

## * Information about the R session used for this document

## chunk 52
sessionInfo()


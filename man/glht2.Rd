% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sCorrect-glht2.R
\name{glht2}
\alias{glht2}
\alias{glht2.lvmfit}
\alias{glht2.lvmfit2}
\alias{glht2.mmm}
\alias{glht.lvmfit2}
\title{General Linear Hypothesis Testing With Small Sample Correction}
\usage{
glht2(object, ...)

\method{glht2}{lvmfit}(
  object,
  linfct,
  rhs = NULL,
  robust = FALSE,
  cluster = NULL,
  ssc = lava.options()$ssc,
  df = lava.options()$df,
  ...
)

\method{glht2}{lvmfit2}(object, linfct, rhs = NULL, robust = FALSE, cluster = NULL, ...)

\method{glht2}{mmm}(object, linfct, rhs = 0, robust = FALSE, cluster = NULL, ...)

\method{glht}{lvmfit2}(model, linfct, rhs = NULL, robust = FALSE, cluster = NULL, ...)
}
\arguments{
\item{object, model}{a \code{lvmfit}, \code{lvmfit2}, or \code{mmm} object.}

\item{...}{[logical] arguments passed to lower level methods.}

\item{linfct}{[matrix or vector of character] the linear hypotheses to be tested. Same as the argument \code{par} of \code{\link{createContrast}}.}

\item{rhs}{[vector] the right hand side of the linear hypotheses to be tested.}

\item{robust}{[logical] should robust standard error be used? 
Otherwise rescale the influence function with the standard error obtained from the information matrix.}

\item{cluster}{[integer vector] the grouping variable relative to which the observations are iid.}

\item{ssc}{[character] method used to correct the small sample bias of the variance coefficients: no correction (\code{"none"}/\code{FALSE}/\code{NA}),
correct the first order bias in the residual variance (\code{"residual"}), or correct the first order bias in the estimated coefficients \code{"cox"}).
Only relevant when using a \code{lvmfit} object.}

\item{df}{[character] method used to estimate the degree of freedoms of the Wald statistic: Satterthwaite \code{"satterthwaite"}. 
Otherwise (\code{"none"}/\code{FALSE}/\code{NA}) the degree of freedoms are set to \code{Inf}.
Only relevant when using a \code{lvmfit} object.}
}
\value{
A \code{glht} object.
}
\description{
Test linear hypotheses on coefficients from a latent variable models with small sample corrections.
}
\details{
Whenever the argument linfct is not a matrix, it is passed to the function \code{createContrast} to generate the contrast matrix and, if not specified, rhs. \cr \cr

Since only one degree of freedom can be specify in a glht object and it must be an integer, the degree of freedom of the denominator of an F test simultaneously testing all hypotheses is retained, after rounding. \cr \cr

Argument rhs and null are equivalent.
This redondance enable compatibility between \code{lava::compare}, \code{compare2}, \code{multcomp::glht}, and \code{glht2}.
}
\examples{
library(multcomp)

## Simulate data
mSim <- lvm(c(Y1,Y2,Y3)~ beta * eta, Z1 ~ E, Z2 ~ E, Age[40:5]~1)
latent(mSim) <- "eta"
set.seed(10)
n <- 1e2

df.data <- lava::sim(mSim, n, latent = FALSE, p = c(beta = 1))

#### Inference on a single model ####
e.lvm <- estimate(lvm(Y1~E), data = df.data)
summary(glht2(e.lvm, linfct = c("Y1~E + Y1","Y1")))

#### Inference on separate models ####
## fit separate models
lvmX <- estimate(lvm(Z1 ~ E), data = df.data)
lvmY <- estimate(lvm(Z2 ~ E + Age), data = df.data)
lvmZ <- estimate(lvm(c(Y1,Y2,Y3) ~ eta, eta ~ E), 
                 data = df.data)

#### create mmm object #### 
e.mmm <- mmm(X = lvmX, Y = lvmY, Z = lvmZ)

#### create contrast matrix ####
resC <- createContrast(e.mmm, linfct = "E")

#### adjust for multiple comparisons ####
e.glht2 <- glht2(e.mmm, linfct = c(X="E"), df = FALSE)
summary(e.glht2)

}
\seealso{
\code{\link{createContrast}} to create contrast matrices. \cr
\code{\link{estimate2}} to pre-compute quantities for the small sample correction.
}
\concept{multiple comparison}
\concept{multiple comparisons}

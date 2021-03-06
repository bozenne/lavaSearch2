% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sCorrect.R
\name{sCorrect}
\alias{sCorrect}
\alias{sCorrect.lm}
\alias{sCorrect.lm2}
\alias{sCorrect.gls}
\alias{sCorrect.gls2}
\alias{sCorrect.lme}
\alias{sCorrect.lme2}
\alias{sCorrect.lvmfit}
\alias{sCorrect.lvmfit2}
\alias{sCorrect<-}
\alias{sCorrect<-.lm}
\alias{sCorrect<-.lm2}
\alias{sCorrect<-.gls}
\alias{sCorrect<-.gls2}
\alias{sCorrect<-.lme}
\alias{sCorrect<-.lme2}
\alias{sCorrect<-.lvmfit}
\alias{sCorrect<-.lvmfit2}
\title{Satterthwaite Correction and Small Sample Correction}
\usage{
sCorrect(object, adjust.Omega, adjust.n, df, numeric.derivative, param,
  data, tol, n.iter, trace, ...)

\method{sCorrect}{lm}(object, adjust.Omega = TRUE, adjust.n = TRUE,
  df = TRUE, numeric.derivative = FALSE, param = NULL, data = NULL,
  tol = 1e-05, n.iter = 20, trace = 0, ...)

\method{sCorrect}{lm2}(object, ...)

\method{sCorrect}{gls}(object, adjust.Omega = TRUE, adjust.n = TRUE,
  df = TRUE, numeric.derivative = FALSE, param = NULL, data = NULL,
  tol = 1e-05, n.iter = 20, trace = 0, cluster, ...)

\method{sCorrect}{gls2}(object, ...)

\method{sCorrect}{lme}(object, adjust.Omega = TRUE, adjust.n = TRUE,
  df = TRUE, numeric.derivative = FALSE, param = NULL, data = NULL,
  tol = 1e-05, n.iter = 20, trace = 0, cluster, ...)

\method{sCorrect}{lme2}(object, ...)

\method{sCorrect}{lvmfit}(object, adjust.Omega = TRUE, adjust.n = TRUE,
  df = TRUE, numeric.derivative = FALSE, param = NULL, data = NULL,
  tol = 1e-05, n.iter = 20, trace = 0, ...)

\method{sCorrect}{lvmfit2}(object, ...)

sCorrect(x, ...) <- value

\method{sCorrect}{lm}(x, ...) <- value

\method{sCorrect}{lm2}(x, ...) <- value

\method{sCorrect}{gls}(x, ...) <- value

\method{sCorrect}{gls2}(x, ...) <- value

\method{sCorrect}{lme}(x, ...) <- value

\method{sCorrect}{lme2}(x, ...) <- value

\method{sCorrect}{lvmfit}(x, ...) <- value

\method{sCorrect}{lvmfit2}(x, ...) <- value
}
\arguments{
\item{object, x}{a \code{gls}, \code{lme}, or \code{lvm} object.}

\item{adjust.Omega}{[logical] should the standard errors of the coefficients be corrected for small sample bias?}

\item{adjust.n}{[logical] should the correction for the degree of freedom be performed?}

\item{df}{[logical] should the degree of freedoms of the Wald statistic be computed using the Satterthwaite correction?
Otherwise the degree of freedoms are set to \code{Inf}, i.e. a normal distribution is used instead of a Student's t distribution when computing the p-values.}

\item{numeric.derivative}{[logical] should a numerical derivative be used to compute the first derivative of the information matrix?
Otherwise an analytic formula is used.}

\item{param}{[numeric vector, optional] the values of the parameters at which to perform the correction.}

\item{data}{[data.frame, optional] the dataset relative to which the correction should be performed.}

\item{tol}{[numeric >0] the minimum absolute difference between two estimation of the small sample bias.
Below this value, the algorithm used to estimate the bias stop.}

\item{n.iter}{[integer >0] the maximum number of iterations used to estimate the small sample bias of the residual variance-covariance matrix.}

\item{trace}{[logical] should the execution of the function be traced.}

\item{...}{[internal] only used by the generic method or by the <- methods.}

\item{cluster}{[integer vector] the grouping variable relative to which the observations are iid.
Only required for \code{gls} models with no correlation argument.}

\item{value}{[logical] value for the arguments \code{adjust.Omega} and \code{adjust.n}.}
}
\description{
Correct the bias of the ML estimate of the variance and compute the first derivative of the information matrix.
}
\details{
The argument \code{value} is equivalent to the argument \code{bias.correct} of the function \code{summary2}.
}
\examples{
n <- 5e1
p <- 3
X.name <- paste0("X",1:p)
link.lvm <- paste0("Y~",X.name)
formula.lvm <- as.formula(paste0("Y~",paste0(X.name,collapse="+")))

m <- lvm(formula.lvm)
distribution(m,~Id) <- sequence.lvm(0)
set.seed(10)
d <- lava::sim(m,n)

## linear model
e.lm <- lm(formula.lvm,data=d)
system.time(
sCorrect(e.lm) <- TRUE ## i.e. bias.correct = TRUE
)

## gls model
library(nlme)
e.gls <- gls(formula.lvm, data = d, method = "ML")
sCorrect(e.gls, cluster = 1:NROW(d)) <- TRUE ## i.e. bias.correct = TRUE
summary2(e.gls)

## latent variable model
e.lvm <- estimate(lvm(formula.lvm),data=d)
sCorrect(e.lvm) <- TRUE ## i.e. bias.correct = TRUE
summary2(e.lvm)

}
\concept{derivative of the score equation}
\concept{small sample inference}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onestepADMM.R
\name{onestepADMM}
\alias{onestepADMM}
\title{ADMM function}
\usage{
onestepADMM(
  x,
  y,
  family = c("gaussian", "cox", "poisson", "logistic"),
  rho = 1,
  censor = NULL,
  C = NULL,
  d = NULL,
  inexact = TRUE,
  initial.guess = NULL,
  tol = 1e-07,
  maxit = 1000,
  lambda = NULL,
  intercept = F,
  min.abs = 0,
  lower.limits = -Inf,
  upper.limits = Inf,
  penalty.factor = 1,
  delta_0 = NULL,
  u1_0 = NULL,
  u2_0 = NULL,
  z_0 = NULL,
  hessian = NULL,
  exclude = rep(FALSE, nvars),
  trace = T
)
}
\arguments{
\item{x}{Design matrix, data for the predictors}

\item{y}{Response vector}

\item{family}{The type of model to fit; at this point,
only logistic and cox models have been updated with the subset approach outlined in the paper}

\item{rho}{proximal mapping parameter for z-update in ADMM}

\item{censor}{If a cox model is fit, the censoring vector for the survival data}

\item{C}{the matrix that linearly constrains the regression coefficients, e.g. a 1 by p
vector constrains the sum of the coefficients}

\item{d}{the value that the constrained regression coefficients should equal e.g. the above
constraint vector and d = 0 would be a sum-to-zero constraint}

\item{inexact}{Whether ADMM should be fit with an inexact speedup, required for all non-gaussian models}

\item{initial.guess}{initial regression coefficient guess}

\item{tol}{tolerance for ADMM, determining primal and dual conditions for termination}

\item{maxit}{maximum ADMM iterations, required do the slowness of ADMM for high accuracy convergence}

\item{lambda}{Lasso penalty parameter}

\item{intercept}{Whether an unpenalized, unconstrained intercept should be included in the model}

\item{min.abs}{The minimum value such that terms smaller in absolute magnitude are shrunk to 0}

\item{lower.limits}{If values are less than lower.limits, they are set to lower.limits}

\item{upper.limits}{If values are greater than upper.limits, they are set to upper.limits}

\item{penalty.factor}{The vector that multiplies the lasso penalty, allowing for different penalties
for each term e.g. 1 is penalization equal to lambda and 0 is no penalization}

\item{delta_0}{output from prior ADMM as new start point for rerun ADMM}

\item{u1_0}{output from prior ADMM as new start point for rerun ADMM}

\item{u2_0}{output from prior ADMM as new start point for rerun ADMM}

\item{z_0}{output from prior ADMM as new start point for rerun ADMM}

\item{hessian}{output from prior ADMM as new start point for rerun ADMM}

\item{exclude}{predictors to exclude from the model, coefficient forced to equal 0}

\item{trace}{if true, also outputs time to fit, solution at each ADMM iteration, Hessian matrix, and
signs of the regression coefficients}
}
\description{
This function allows you to use ADMM to fit the constrained lasso
}
\examples{
onestepADMM()
}
\keyword{ADMM}

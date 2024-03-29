% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inner_solution.R
\name{inner.solution}
\alias{inner.solution}
\title{Function for constrained solution of subset (inner solution)}
\usage{
inner.solution(
  admm.fit,
  family = c("logistic", "cox"),
  x,
  y = NULL,
  lambda,
  censor = NULL,
  tolerance = 1e-07,
  intercept = F,
  C,
  penalty.factor,
  d = 0,
  trace = F
)
}
\arguments{
\item{admm.fit}{The ADMM fit outputted by onestepADMM}

\item{family}{The type of model. Note for now only logistic and Cox models are supported.}

\item{x}{The design matrix}

\item{y}{The response vector}

\item{lambda}{The lasso penalty}

\item{censor}{If using a Cox model, the censoring vector}

\item{tolerance}{Tolerance of ADMM, influences the primal and dual criteria for termination}

\item{intercept}{Whether an unpenalized an unconstrained intercept should be added to the model}

\item{C}{The constraint that linearly combines the regression coefficients. Should be a matrix}

\item{penalty.factor}{The penalization vector that multiplies the lasso penalty lambda e.g. elements = 1
are lasso penalized and elements = 0 are unpenalized}

\item{d}{The scalar value that C times the regrssion coefficients should equal}

\item{trace}{If false, returns regression coefficients. If true, also returns number of ADMM iterations required}
}
\description{
This function uses repeated ADMM and Newton-Raphson to find a solution for the predictor subset
}
\keyword{inner_solution}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KKTtestfct.R
\name{KKTtestfct}
\alias{KKTtestfct}
\title{Check if KKT conditions hold}
\usage{
KKTtestfct(
  par_guess,
  family = c("logistic", "cox"),
  lambda,
  x,
  y = NULL,
  censor = NULL,
  C = matrix(1, nrow = 1, ncol = ncol(x))
)
}
\arguments{
\item{par_guess}{The guess of the regression coefficients}

\item{family}{The type of model. For now, only logistic and Cox models are supported}

\item{lambda}{The value of the lasso penalty}

\item{x}{The design matrix}

\item{y}{The response vector}

\item{censor}{If using a Cox model, the censoring vector}

\item{C}{The linear constraint, as a 1 by p matrix}

\item{trace}{If false, returns true/false output for KKT conditions. If true, returns that and the following:
whether the KKT conditions hold for each of unconstrained and constrained terms, and the eta vectors.
Note for KKT to hold the max of the eta_lower vector must be leq the min of the eta_upper vector}
}
\description{
This function checks if the proposed solution satisfies the KKT conditions for optimality
}
\keyword{KKTtest}

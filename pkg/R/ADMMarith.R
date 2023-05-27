#' Basic math functions for ADMM
#' 
#' @param x numeric vector
#' @param kappa shrinkage factor
#' @name ADMMarith
NULL
#> NULL

#' @rdname ADMMarith
vector.2.norm <- function(x) {
  sqrt(sum(x ^ 2))
}

#' @rdname ADMMarith
shrinkage <- function(x, kappa){
  
  z = pmax(0, x - kappa) - pmax(0, -x - kappa)
  
}

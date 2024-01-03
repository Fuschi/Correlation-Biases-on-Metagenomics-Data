#' Get upper triangular values
#' 
#' @description Returns as vector the upper triangular values of a square symmetric 
#' matrix
#' 
#' @param X numeric symmetric matrix 
#' 
TRIU <- function(X){
  
  if(!is.matrix(X)) stop("X must be a matrix or a vector")
  if(!is.numeric(X)) stop("X must be numeric")
  if(!isSymmetric(X)) stop("X must be symmetric")
  
  return(X[upper.tri(X)])
}
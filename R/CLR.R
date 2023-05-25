#' Centered Log-Ratio Transformation
#' 
#' @description It calculates the centered log-ratio transformation of X. 
#' The zeroes are replaced with 0.65 (https://doi.org/10.1016/j.chemolab.2021.104248).
#'
#' @param X numeric matrix with all elements greater than or equal to 0.
#' @param mar Integer giving the dimension where the function will be applied;
#' 1 for rows and 2 for columns (default 1).
#' 
#' @export
CLR <- function(X, mar=1){
  
  if(!is.matrix(X)) stop("X must be a matrix or a vector")
  if(!is.numeric(X) | any(X<0)) stop("X must be numeric with all elements greater than or equal to 0")
  if(!(mar%in%c(1,2))) stop("mar has as possible values only 1 and 2.")
  
  X[X==0] <- .65
  ref <- apply(X, mar, function(x) mean(log(x)) )
  return(as.matrix(log(X) - ref))
  
}
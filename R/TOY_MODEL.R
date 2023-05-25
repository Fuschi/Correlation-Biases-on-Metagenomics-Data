#' Meta-genomic Toy Model
#' 
#' @description numerical simulations for meta-genomic data with user-chosen 
#' correlation structure and sample heterogeneity.
#'
#' @param n number of observations.
#' @param cor correlation structure of resulting data; cor must be a number or a
#' correlation matrix. If it is a number it is understood that it is the 
#' correlation between two variables.
#' @param D alternative choice to cor and indicates dimensionality. Produces a 
#' correlation matrix of dimension D with by default two pairs of strongly 
#' correlated variables, one positively and one negatively with correlations 
#' values equal to +0.8 and -0.8. D cannot take values less than 5.
#' @param M magnification factor, real positive number that modify the 
#' heterogeneity of the samples. In practice the function multiplies by a 
#' variable by the value of M.
#' @param dist characters indicate the chosen distribution.
#' @param param list containing the parameters of the distribution for each 
#' dimension. All elements in the list must be named as the parameter name of
#' the distribution. If the elements have all lengths 1 the values are repeated
#' for all the dimensions.
#' @param method type of correlation used. Possible choices are "pearson", 
#' "kendall", "spearman" (default pearson). 
#' @param seed random seed for reproducibility (default 42).
#' @param force.positive logical indicates when to force all elements of Y
#' to be positive number adding the minimum to all others (this passage does not
#' affect correlations).
#'
#' @returns \code{toy_model} returns and object of class "toy_model" containing:
#' \itemize{
#'  \item cor input correlation matrix.
#'  \item normal generated gaussian data.
#'  \item cor_normal correlation matrix of normal.
#'  \item NorTA generated data of choosen distribution.
#'  \item cor_NorTA correlation matrix of Y.
#'  \item L1 relative abundances of NorTA.
#'  \item cor_L1 correlation matrix of L1.
#'  \item CLR matrix of abundances of clr transformed data from NorTA.
#'  \item cor_CLR correlation matrix of CLR.
#' }
#' @export
TOY_MODEL <- function(n,cor,D,M,dist,param=list(),method="pearson",seed=42,
                      force.positive=FALSE){
  
  #check n
  #-----------------------------------#
  if(!is.numeric(n) | n <=0 | round(n)!=n) stop("n must be positve integer number")
  #check dist
  #-----------------------------------#
  ifelse(is.character(dist), 
         qdist<-paste("q",dist,sep=""),
         stop("dist must be a character string"))
  if(!exists(qdist))stop(paste("quantile function",qdist,
                               "has not been found, maybe you need to include the required package?"))
  #check cor and D
  #-----------------------------------#
  if(missing(cor) && missing(D)) stop("cor and D cannot be unassigned togheter")
  if(!missing(cor) && !missing(D)) stop("cor and D cannot be assigned togheter")
  #check cor
  #-----------------------------------#
  if(!missing(cor)){
    if(!is.numeric(cor) | !isSymmetric(cor) | any(diag(cor)!=1) | any(abs(cor)>1))
      stop("cor must be a symmetric matrix with all elements in range [-1,1] and 
           with all elements values equal to 1")
  }
  #check D
  #-----------------------------------#
  if(!missing(D)){
    if(!is.numeric(D) | D<5 | round(D)!=D) stop("D must be positve integer number greater or equal to 5")
  }
  #check param
  #-----------------------------------#
  if(length(param)!=0){
    lengths <- lapply(param, length)
    if(!all(lengths==1 | lengths==D)){
      stop("all elemnts in param must have lengts equal to 1 or equal to the dimension of cor")
    } 
    param <- lapply(param, function(x) if(length(x)==1) rep(x,D) else x )
  }
  #check M
  #-----------------------------------#
  if(!is.numeric(M) | M<=0 | round(M)!=M) stop("M must be positve number")
  #check method
  #-----------------------------------#
  method  <- match.arg(method,c("pearson", "kendall", "spearman"))
  
  
  if(missing(cor)){
    cor <- matrix(0,nrow=D,ncol=D)
    cor[ceiling(.21*D),ceiling(.61*D)] <- .8
    cor[ceiling(.41*D),ceiling(.81*D)] <- -.8
    cor <- cor+t(cor)
    diag(cor) <- 1
  } else {
    D <- nrow(cor)
  }
  
  # Generate random Gaussian variables with custom correlation structure
  normal <- mvtnorm::rmvnorm(n=n,sigma=cor)
  # Get probabilities
  unif <- stats::pnorm(normal)
  
  # Transform data to wanted distribution 
  NorTA <- sapply(1:nrow(cor), function(idx){
    sub.param <- c(list(p=unif[,idx]),lapply(param, `[[`, idx))
    do.call(what=paste("q",dist,sep=""),args=sub.param)
  })
  
  # Move all elements to positive values if required
  if(force.positive && any(NorTA<0)) NorTA <- NorTA - min(NorTA)
  
  # Apply the magnification factor to the first variable
  NorTA[,1] <- NorTA[,1]*M
  
  # Check that no elements are negative
  if(any(NorTA<0)) stop("the transformed data NorTA cannot have negative values.")
  
  # Get correlation of generated data and transformed NorTA data
  cor_normal = stats::cor(normal,method=method)
  cor_NorTA = stats::cor(NorTA,method=method)
  
  # Make L1 normalization and calculate correlation
  L1 <- NorTA/rowSums(NorTA)
  cor_L1 <- stats::cor(L1,method=method)
  
  # clr normalization
  clr <- function(X){
    if(any(X==0)) X <- X+1
    ref <- apply(X, 1, function(x) mean(log(x)) )
    return(as.matrix(log(X) - ref))
  }
  
  # Make CLR normalization and calculate correlation
  CLR <- clr(NorTA)
  cor_CLR <- stats::cor(CLR,method=method)
  
  # Resume results in an object
  results <- list()
  results$cor <- cor
  results$normal <- normal
  results$cor_normal <- cor_normal 
  results$NorTA <- NorTA
  results$cor_NorTA <- cor_NorTA
  results$L1 <- L1
  results$cor_L1 <- cor_L1
  results$CLR <- CLR
  results$cor_CLR <- cor_CLR
  class(results) <- "toy_model"
  
  return(results)
}
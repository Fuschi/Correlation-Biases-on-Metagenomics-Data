#' Hurdle Truncated Log-Normal Model
#' 
#' @description mle parameter estimation, random generation, density and 
#' quantile function for the hurdle truncated log-normal distribution.
#'
#' @param x vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param phi probability of zero.
#' @param meanlog mean in the log scale.
#' @param sdlog standard deviation in log scale.
#' @param a lower bounds.
#' @param b upper bounds.
#' @param y numeric vector containing only finite values.
#' @param integer (default TRUE) logical, indicates when to approx the results 
#' to integer.
#' @param warning.silent logical; if TRUE suppress all warning messages from 
#' \code{\link{fitdistr}} during the mle process.


mle.htrlnorm <- function(y, warning.silent=TRUE){
  
  if(!is.vector(y) | any(y<0) | !is.numeric(y)) stop("y must be a numeric vector with all values >= 0.")
  
  y <- log1p(y)
  y1 <- y[y>0]
  
  n <- length(y)
  n1 <- length(y1)
  n0 <- n - n1
  
  phi <- n0/n
  
  param.trnorm.mle <- NULL
  try(param.trnorm.mle <- MASS::fitdistr(x=y1, densfun=truncnorm::dtruncnorm,
                                         start=list(mean=mean(y1),sd=stats::sd(y1)),
                                         a=0,b=Inf),
      silent=warning.silent
  )
  
  if(!is.null(param.trnorm.mle) && phi!=0){
    loglik <- n0*log(phi) + n1*log(1-phi) + param.trnorm.mle$loglik
  }else if (!is.null(param.trnorm.mle) && phi==0){
    loglik <- param.trnorm.mle$loglik
  }else{
    loglik <- NA
  }
  
  
  if(!is.null(param.trnorm.mle)){
    estimate <- c("phi"=phi,"meanlog"=as.numeric(param.trnorm.mle$estimate[1]),
                  "sdlog"  =as.numeric(param.trnorm.mle$estimate[2]))
    result <- list("estimate"=estimate,
                   "loglik" =loglik,
                   "success"= 1)
  } else {
    estimate <- c("phi"=phi,"meanlog"=NA,"sdlog"=NA)
    result <- list("estimate"=estimate,
                   "loglik" =loglik,
                   "success"= 1)
  }
  
  return(result)
  
}

dhtrlnorm <- function(x, phi=0, meanlog=0, sdlog=1, a=0, b=Inf){
  
  #CHECK ARGUMENTS
  if(!is.vector(x) | !is.numeric(x)) stop("x must be a numeric vector")
  if(!is.numeric(phi) | !is.numeric(meanlog) | !is.numeric(sdlog)) stop("phi, meanlog, sdlog must be all number")
  if(phi<0 || phi>1) stop("phi must be in range [0,1]")
  if(sdlog<=0) stop("sdlog must be greater than 0")
  if(!is.numeric(a)) stop("a must be numeric")
  if(!is.numeric(b)) stop("b must be numeric")
  
  ans <- (1-phi)*dtruncnorm(x,a=a,b=b,mean=meanlog,sd=sdlog)
  ans[x==0] <- phi
  
  return(ans)
}

qhtrlnorm <- function(p, phi=0, meanlog=0, sdlog=1, a=0, b=Inf, integer=TRUE){
  
  #CHECK ARGUMENTS
  if(!is.vector(p) | !is.numeric(p)) stop("p must be a numeric vector")
  if(!is.numeric(phi) | !is.numeric(meanlog) | !is.numeric(sdlog)) stop("phi, meanlog, sdlog must be all number")
  if(phi<0 || phi>1) stop("phi must be in range [0,1]")
  if(sdlog<0) stop("sdlog must be greater than 0")
  if(!is.numeric(a)) stop("a must be numeric")
  if(!is.numeric(b)) stop("b must be numeric")
  if(!is.logical(integer)) stop("integer must be logical")
  
  if(phi==1){
    ans <- rep(0,length(p))
  } else {
    ans <- rep(NA_real_,length(p))
    ans[p<=phi] <- 0
    
    idx <- is.na(ans)
    ans[idx]  <- qtruncnorm(p=(p[idx]-phi)/(1-phi), a=a, b=b,
                            mean=meanlog, sd=sdlog)
    
    ans <- expm1(ans)
    
    if(integer){
      ans[intersect(which(ans>0),which(ans<1))] <- 1
      ans <- round(ans)
    }
  }
  return(ans)
}

rhtrlnorm <- function(n, phi=0, meanlog=0, sdlog=1, a=0, b=Inf, integer=TRUE){
  
  if(!is.numeric(n) | !is.numeric(phi) | !is.numeric(meanlog) | !is.numeric(sdlog)) stop("n, phi, meanlog, sdlog must be all number")
  if(round(n)!=n) stop("n must be an integer")
  if(n<1) stop("n must be a positive integer number")
  if(phi<0 || phi>1) stop("phi must be a number in [0,1]")
  if(sdlog<0) stop("sdlog must be a positive number")
  if(!is.numeric(a)) stop("a must be numeric")
  if(!is.numeric(b)) stop("b must be numeric")
  if(!is.logical(integer)) stop("integer must be logical")
  
  ans <- stats::rbinom(n=n, size=1, prob=1-phi)
  m <- length(ans[ans>0])
  ans[ans==1] <- truncnorm::rtruncnorm(n=m, a=a, b=b,
                                       mean=meanlog,sd=sdlog)
  
  ans <- expm1(ans)
  
  if(integer){
    ans[intersect(which(ans>0),which(ans<1))] <- 1
    ans <- round(ans)
  }
  
  return(ans)
}
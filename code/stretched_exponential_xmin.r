dstexpo_xmin <- function(x, beta = 1, lambda = 1, xmin = 0, log = FALSE){
  if(x <= xmin){
    dens <- -Inf
  }else{
    lc <-  log(beta) + log(lambda) + lambda * xmin^beta
    dens <- lc  + (beta-1)*log(x) - lambda * x^beta
  }
  if(!log) dens <- exp(dens)
  return(dens)
}
dstexpo_xmin <- Vectorize(dstexpo_xmin)
###
pstexpo_xmin <- function(x, beta = 1, lambda = 1, xmin = 0){
  if(x <= xmin){
    ans <- 0
  }else{
    int <- integrate(function(y) dstexpo_xmin(x = y, beta = beta, lambda = lambda, xmin = xmin), xmin, x)
    ans <- int$value 
  }
  return(ans)
}
pstexpo_xmin <- Vectorize(pstexpo_xmin)
##
pstexp <- function(x, beta = 1, lambda = 1, logp = FALSE){
  (1- exp(-lambda * x^beta))
}
##
rstexpo <- function(n, beta, lambda, xmin){
  # TODO
}
##
stexpo_xmin_cdf <- function(x, beta, lambda, xmin){
  if(x < xmin){
    ans <- 0
  }else{
    Fm <- pstexp(x = xmin, beta, lambda)
    Fx <- pstexp(x,beta, lambda) 
    ans <- 1-(Fx-Fm)/(1-Fm)
  }
  return(ans)
} 
stexpo_xmin_cdf <- Vectorize(stexpo_xmin_cdf)
##
library(reticulate)
log_sum_exp <- py_func(logSumExp) ## depends on bplaw_aux.r
log_diff_exp <- py_func(logDiffExp)
##
get_marginal_likelihood_str_expo <- function(x, a1 = 1, b1 = 1, a2 = 2, b2 = 2, xmin = 0){
  N <- length(x)
  sumLog <- sum(log(x))
  lconst_par <- ( a1*log(b1) - lgamma(a1) ) + ( a2*log(b2) - lgamma(a2) )
  lconst_p <- -sumLog
  lconst_N <- lgamma(N + a2)
  fb <- function(bb, args){
    dummy <- args[1] + args[2]
    lSb <-  log_sum_exp(bb * log(x))
    ld <- (N + a1 - 1)*log(bb) - (b1 - sumLog)*bb - (N + a2) * log_diff_exp(log_sum_exp(c(lSb, log(b2))), log(N) + bb * log(xmin))
    return(dummy * ld)
  }
  fb <- Vectorize(fb)
  require(reticulate)
  lint <- import("lintegrate", convert = FALSE)
  intt <- lint$lcquad(py_func(fb), a = r_to_py(0), b = r_to_py(5), ## increasing this bound will actually DECREASE accuracy
                      args=c(1, 0), epsabs = r_to_py(1e-30))
  res <- reticulate::py_to_r(intt[0])
  return(list(logml  = lconst_par + lconst_N + lconst_p + res, 
              otherstuff = intt ))
}
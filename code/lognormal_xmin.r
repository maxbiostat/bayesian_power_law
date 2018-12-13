erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
dlnorm_xmin <- function(x, meanlog = 0, sdlog = 1, xmin = 0, log = FALSE){
  if(x < xmin){
    dens <- -Inf
  }else{
    lc <-  log(2) - log(erfc( (log(xmin)-meanlog)/(sqrt(2)*sdlog) ))
    dens <- lc + dlnorm(x, meanlog, sdlog, log = TRUE) 
  }
  if(!log) dens <- exp(dens)
  return(dens)
}
dlnorm_xmin <- Vectorize(dlnorm_xmin)
##
plnorm_xmin <- function(x, meanlog = 0, sdlog = 1, xmin = 0, log = FALSE){
  if(x < xmin){
    ans <- 0
  }else{
   intt <- integrate(function(y) dlnorm_xmin(x = y, meanlog = meanlog, sdlog = sdlog, xmin = xmin), 0, x)
   ans <- intt$value
  }
  return(ans)
}
plnorm_xmin <- Vectorize(plnorm_xmin)
#####################################
#####################################
# m <- .4
# mu <- -.8
# sigma <- .3
# integrate(function(x) dlnorm_xmin(x, meanlog = mu, sdlog = sigma, xmin = m), 0, Inf)
# curve(dlnorm_xmin(x, meanlog = mu, sdlog = sigma, xmin = m))
# #####
# kappa <- function(mu, sigma, xm){
#   sqrt(2/(pi * sigma^2)) *  1/erfc( (log(xm)-mu)/(sqrt(2)*sigma) )
# }
# #####
# X <- rlnorm(1000, mu, sigma)
# Y <- X[X > m]
# samp <- Y[1:10]
# n <- length(samp)
# Sp <- sum(log(samp))
# 
# sum(
#   dlnorm_xmin(samp, mu, sigma, m, log = TRUE)
# )
# n*log(kappa(mu, sigma, m)) - (n *mu^2)/(2*sigma^2) + (mu/sigma^2 -1)*Sp - sum(log(samp)^2)/(2*sigma^2)
compress_data <- function(x){
  tab <- table(x)
  return(data.frame(x = as.numeric(names(tab)), n = as.numeric(tab)))
}
##
logSumExp <- function(x) log(sum(exp(x - max(x)))) + max(x)
##
logDiffExp <- function(x, y){
  biggest <- x
  xx1 = x - biggest
  xx2 = y - biggest
  ans <- log(exp(xx1) - exp(xx2)) + biggest
  return(ans)
}
##
lgamma_inc <- function(s, x) pgamma(x, s, lower=FALSE, log.p = TRUE) + lgamma(s)
# OLS: returns cbind(coef, se) 
ols <- function(y, x) {
  
  b <- solve(crossprod(x), crossprod(x, y))
  e <- y - x %*% b
  s2 <- as.numeric(crossprod(e)) / (nrow(x) - ncol(x))
  se <- sqrt(diag(s2 * solve(crossprod(x))))
  
  cbind(b, se)
  
}

# R2OS statistic + Clark-West (2007) one-sided p-value 
# Returns c(R2OS, p-value) where R2OS is in percent.
# H0: MSFE_0 <= MSFE_i  vs  H1: MSFE_0 > MSFE_i (i.e. R2OS > 0)
r2os_stat <- function(e_ha, e_alt, fc_ha, fc_alt) {
  
  r2 <- 100 * (1 - sum(e_alt ^ 2) / sum(e_ha ^ 2))
  adj <- e_ha ^ 2 - (e_alt ^ 2 - (fc_ha - fc_alt) ^ 2)
  tst <- coef(summary(lm(adj ~ 1)))[, "t value"]
  
  c(r2, 1 - pnorm(tst))
  
}

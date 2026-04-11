# Returns list(util, w): mean-variance utility and risky weights.
# Weights are clamped to [0, 1.5].
asset_alloc <- function(y, rf, fc, vol, gamma) {
  
  w   <- pmax(0, pmin(1.5, (1 / gamma) * (fc / vol)))
  
  ret <- rf + w * y
  
  list(util = mean(ret) - 0.5 * gamma * var(ret), w = w)
  
}

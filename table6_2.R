#'---
#' title: Table 6.2 – Monthly U.S. Equity Premium Out-of-Sample Forecasting
#'         Results Based on Multiple Economic Variables, 1957:01–2010:12
#'
#' Rapach & Zhou (2013), Chapter 6, Handbook of Economic Forecasting Vol. 2A
#'
#' Run refs/data_prep.R first to generate data/refs_datasets.rds.
#'
#' Panel A: Unrestricted forecasts
#' Panel B: Forecasts with non-negativity restrictions
#' Methods: Kitchen sink | SIC | POOL-AVG | POOL-DMSFE | Diffusion index |
#'          Sum-of-the-parts
#'
#' R2OS computed with LOG equity premium (paper p.353).
#' Utility gain (Delta) computed with SIMPLE equity premium (footnote 28).
#'---

# Read packages
library("tidyverse")

# Read auxiliary functions
source("R/asset-alloc.R")
source("R/clark-west.R")

ds        <- readRDS("data/rz2013_data.rds")
data_log  <- ds$data_log    # log eq premium + 14 predictors
data_sim  <- ds$data_sim    # simple eq premium + 14 predictors
sink_data <- ds$sink_data   # 12-predictor kitchen-sink
sop_data  <- ds$sop_data    # SOP variables (log & sim versions)
rf_lag    <- ds$rf_lag
rec       <- ds$rec

# Setup
y_log <- data_log$eq_prem
y_sim <- data_sim$eq_prem
x     <- select(data_log, -eq_prem)    # 14 individual predictors
xs    <- select(sink_data, -eq_prem)   # 12 kitchen-sink predictors
n     <- nrow(data_log)
nc    <- ncol(x)    # 14
ncs   <- ncol(xs)   # 12

# Pre-convert to matrix for speed
x_mat  <- as.matrix(x)
xs_mat <- as.matrix(xs)

# Full-sample beta signs (log premium)
bf <- vapply(seq_len(nc), function(j) {
  
  ols(y_log[2:n], cbind(x_mat[1:(n - 1), j], 1))[1, 1]
  
}, numeric(1))

bf[5] <- 1  # SVAR sign override

r_win  <- (1946 - 1926) * 12 + 1
p0     <- (1956 - 1946) * 12
p      <- n - r_win - p0
np     <- p0 + p
ma_sop <- 20 * 12 # 20-year MA for SOP earnings growth
theta  <- 0.75 # POOL-DMSFE discount factor
r_di   <- 1L # number of principal components for DI

# Pre-compute SIC model index sets: all subsets of size 1, 2, 3 from 12 sink predictors
sic_combos <- c(
  combn(ncs, 1, simplify = FALSE),
  combn(ncs, 2, simplify = FALSE),
  combn(ncs, 3, simplify = FALSE)
)

n_sic <- length(sic_combos) # 12 + 66 + 220 = 298

# Forecast matrices
fc_ha_log <- numeric(np)
fc_ec_log <- matrix(0, np, nc) # individual log-premium (for POOL-AVG / DMSFE)
fc_ha_sim <- numeric(np)
fc_ec_sim <- matrix(0, np, nc) # individual simple-premium

# Combination forecasts – log premium (for R2OS)
fc_ks_log    <- numeric(np)
fc_sic_log   <- numeric(np)
fc_avg_log   <- numeric(np)
fc_dmsfe_log <- numeric(np)
fc_di_log    <- numeric(np)
fc_sop_log   <- numeric(np) # Forecasts_monthly_log.m formula

# Combination forecasts – simple premium (for Delta)
fc_ks_sim    <- numeric(np)
fc_sic_sim   <- numeric(np)
fc_avg_sim   <- numeric(np)
fc_dmsfe_sim <- numeric(np)
fc_di_sim    <- numeric(np)
fc_sop_sim   <- numeric(np) # Forecasts_monthly.m formula

# Forecast loop
for (i in seq_len(np)) {
  
  ri <- r_win + i - 1
  
  fc_ha_log[i] <- mean(y_log[1:ri])
  fc_ha_sim[i] <- mean(y_sim[1:ri])

  xt  <- x_mat[1:(ri - 1), , drop = FALSE] # (ri-1) x 14
  xts <- xs_mat[1:(ri - 1), , drop = FALSE] # (ri-1) x 12
  yl  <- y_log[2:ri] # (ri-1) log-premium outcomes
  ys  <- y_sim[2:ri] # (ri-1) simple-premium outcomes
  Tk  <- ri - 1L # training sample size

  # Individual predictive regression forecasts (needed for POOL-AVG / DMSFE)
  for (j in seq_len(nc)) {
    
    xj <- cbind(xt[, j], 1)
    xp <- c(x_mat[ri, j], 1)
    fc_ec_log[i, j] <- xp %*% ols(yl, xj)[, 1]
    fc_ec_sim[i, j] <- xp %*% ols(ys, xj)[, 1]
    
  }

  if (i > p0) {

    # Kitchen sink: OLS on all 12 sink predictors
    xts1    <- cbind(xts, 1)
    xp_sink <- c(xs_mat[ri, ], 1)
    fc_ks_log[i] <- xp_sink %*% ols(yl, xts1)[, 1]
    fc_ks_sim[i] <- xp_sink %*% ols(ys, xts1)[, 1]

    # SIC: BIC model selection over 1-/2-/3-predictor subsets
    # Criterion: log(RSS/T) + log(T)*k/T  (k = predictors + intercept)
    sic_vals_log <- numeric(n_sic)
    sic_vals_sim <- numeric(n_sic)
    sic_fc_log <- numeric(n_sic)
    sic_fc_sim <- numeric(n_sic)

    for (km in seq_len(n_sic)) {
      
      sel  <- sic_combos[[km]]
      Xjk  <- cbind(xts[, sel, drop = FALSE], 1)
      xpjk <- c(xs_mat[ri, sel], 1)
      kk   <- length(sel) + 1L
      qr_jk <- qr(Xjk)

      e_l <- qr.resid(qr_jk, yl)
      sic_vals_log[km] <- log(sum(e_l ^ 2) / Tk) + log(Tk) * kk / Tk
      sic_fc_log[km] <- xpjk %*% qr.coef(qr_jk, yl)

      e_s <- qr.resid(qr_jk, ys)
      sic_vals_sim[km] <- log(sum(e_s ^ 2) / Tk) + log(Tk) * kk / Tk
      sic_fc_sim[km] <- xpjk %*% qr.coef(qr_jk, ys)
      
    }
    
    fc_sic_log[i] <- sic_fc_log[which.min(sic_vals_log)]
    fc_sic_sim[i] <- sic_fc_sim[which.min(sic_vals_sim)]

    # POOL-AVG: simple mean of 14 individual forecasts
    fc_avg_log[i] <- mean(fc_ec_log[i, ])
    fc_avg_sim[i] <- mean(fc_ec_sim[i, ])

    # POOL-DMSFE: discounted MSFE-weighted average
    # omega_j ∝ 1 / sum_{s=1}^{i-1} theta^{i-1-s} * e_{s,j}^2
    # Errors accumulated over holdout + past evaluation (steps 1 to i-1)
    past_l <- y_log[(r_win + 1L):(r_win + i - 1L)]
    past_s <- y_sim[(r_win + 1L):(r_win + i - 1L)]
    err_l  <- matrix(past_l, i - 1L, nc) - fc_ec_log[1:(i - 1L), ]
    err_s  <- matrix(past_s, i - 1L, nc) - fc_ec_sim[1:(i - 1L), ]
    w_d <- theta^((i - 2L):0L) # (i-1) weights; most recent = 1
    m_l <- colSums(w_d * err_l ^ 2)
    m_s <- colSums(w_d * err_s ^ 2)
    fc_dmsfe_log[i] <- sum(fc_ec_log[i, ] * (1 / m_l) / sum(1 / m_l))
    fc_dmsfe_sim[i] <- sum(fc_ec_sim[i, ] * (1 / m_s) / sum(1 / m_s))

    # Diffusion index: first principal component of 14 predictors
    # Standardise all ri observations → PCA → regress on first r_di PCs
    x_std   <- scale(x_mat[1:ri, ])
    pca_t   <- prcomp(x_std, center = FALSE, scale. = FALSE)
    F_all   <- pca_t$x # ri x nc scores
    F_train <- F_all[1:(ri - 1L), 1:r_di, drop = FALSE] # training scores
    F_pred  <- F_all[ri, 1:r_di] # prediction score
    cbF     <- cbind(F_train, 1)
    qr_di   <- qr(cbF)
    fc_di_log[i] <- c(F_pred, 1) %*% qr.coef(qr_di, yl)
    fc_di_sim[i] <- c(F_pred, 1) %*% qr.coef(qr_di, ys)

    # Sum-of-the-parts
    # Log version (Forecasts_monthly_log.m): 20yr-MA log earnings growth
    #   + log(1+D/P) - log(1+rf)   → used for R2OS
    fc_sop_log[i] <- mean(sop_data$e_growth_log[(ri - ma_sop + 1):ri]) + sop_data$dp_sop_log[ri] - sop_data$rf_sop_log[ri]

    # Simple version (Forecasts_monthly.m): 20yr-MA simple earnings growth
    #   + D/P - rf   → used for Delta
    fc_sop_sim[i] <- mean(sop_data$e_growth_sim[(ri - ma_sop + 1):ri]) + sop_data$dp_sop_sim[ri] - sop_data$rf_sop_sim[ri]
    
  }
  
}

# Non-negativity constrained versions
fc_ks_log_ct    <- pmax(0, fc_ks_log)
fc_sic_log_ct   <- pmax(0, fc_sic_log)
fc_avg_log_ct   <- pmax(0, fc_avg_log)
fc_dmsfe_log_ct <- pmax(0, fc_dmsfe_log)
fc_di_log_ct    <- pmax(0, fc_di_log)
fc_sop_log_ct   <- pmax(0, fc_sop_log)
fc_ks_sim_ct    <- pmax(0, fc_ks_sim)
fc_sic_sim_ct   <- pmax(0, fc_sic_sim)
fc_avg_sim_ct   <- pmax(0, fc_avg_sim)
fc_dmsfe_sim_ct <- pmax(0, fc_dmsfe_sim)
fc_di_sim_ct    <- pmax(0, fc_di_sim)
fc_sop_sim_ct   <- pmax(0, fc_sop_sim)

# Trim to evaluation period (1957:01 – 2010:12)
idx_ev     <- seq(p0 + 1, np)
actual_log <- y_log[(r_win + p0 + 1):n]
actual_sim <- y_sim[(r_win + p0 + 1):n]
fc_ha_log  <- fc_ha_log[idx_ev]
fc_ha_sim  <- fc_ha_sim[idx_ev]

tr <- function(v) v[idx_ev]

fc_ks_log    <- tr(fc_ks_log);    fc_ks_log_ct    <- tr(fc_ks_log_ct)
fc_sic_log   <- tr(fc_sic_log);   fc_sic_log_ct   <- tr(fc_sic_log_ct)
fc_avg_log   <- tr(fc_avg_log);   fc_avg_log_ct   <- tr(fc_avg_log_ct)
fc_dmsfe_log <- tr(fc_dmsfe_log); fc_dmsfe_log_ct <- tr(fc_dmsfe_log_ct)
fc_di_log    <- tr(fc_di_log);    fc_di_log_ct    <- tr(fc_di_log_ct)
fc_sop_log   <- tr(fc_sop_log);   fc_sop_log_ct   <- tr(fc_sop_log_ct)
fc_ks_sim    <- tr(fc_ks_sim);    fc_ks_sim_ct    <- tr(fc_ks_sim_ct)
fc_sic_sim   <- tr(fc_sic_sim);   fc_sic_sim_ct   <- tr(fc_sic_sim_ct)
fc_avg_sim   <- tr(fc_avg_sim);   fc_avg_sim_ct   <- tr(fc_avg_sim_ct)
fc_dmsfe_sim <- tr(fc_dmsfe_sim); fc_dmsfe_sim_ct <- tr(fc_dmsfe_sim_ct)
fc_di_sim    <- tr(fc_di_sim);    fc_di_sim_ct    <- tr(fc_di_sim_ct)
fc_sop_sim   <- tr(fc_sop_sim);   fc_sop_sim_ct   <- tr(fc_sop_sim_ct)

idx_rec <- which(rec$recession == 1)
idx_exp <- which(rec$recession == 0)

# R2OS statistics (log premium)
e_ha <- actual_log - fc_ha_log

# col order: KS | SIC | POOL-AVG | POOL-DMSFE | DI | SOP
fc_mat_a <- cbind(fc_ks_log, fc_sic_log, fc_avg_log, fc_dmsfe_log, fc_di_log, fc_sop_log)
fc_mat_b <- cbind(fc_ks_log_ct, fc_sic_log_ct, fc_avg_log_ct, fc_dmsfe_log_ct, fc_di_log_ct, fc_sop_log_ct)

nm <- ncol(fc_mat_a) # 6

r2_a <- matrix(0, nm, 6)
r2_b <- matrix(0, nm, 6)

for (j in seq_len(nm)) {
  
  e_a <- actual_log - fc_mat_a[, j]
  e_b <- actual_log - fc_mat_b[, j]

  # Overall
  r2_a[j, 1:2] <- r2os_stat(e_ha, e_a, fc_ha_log, fc_mat_a[, j])
  r2_b[j, 1:2] <- r2os_stat(e_ha, e_b, fc_ha_log, fc_mat_b[, j])
  
  # Expansion
  r2_a[j, 3:4] <- r2os_stat(e_ha[idx_exp], e_a[idx_exp], fc_ha_log[idx_exp], fc_mat_a[idx_exp, j])
  r2_b[j, 3:4] <- r2os_stat(e_ha[idx_exp], e_b[idx_exp], fc_ha_log[idx_exp], fc_mat_b[idx_exp, j])
  
  # Recession
  r2_a[j, 5:6] <- r2os_stat(e_ha[idx_rec], e_a[idx_rec], fc_ha_log[idx_rec], fc_mat_a[idx_rec, j])
  r2_b[j, 5:6] <- r2os_stat(e_ha[idx_rec], e_b[idx_rec], fc_ha_log[idx_rec], fc_mat_b[idx_rec, j])
  
}

# Asset allocation / utility gain (simple premium)
gamma   <- 5
win_vol <- 12 * 5
fc_vol  <- numeric(p)

for (t in seq_len(p)) {
  
  idx_v     <- (r_win + p0 + t - win_vol):(r_win + p0 + t - 1)
  yv        <- y_sim[idx_v]
  fc_vol[t] <- mean(yv ^ 2) - mean(yv) ^ 2
  
}

rf_p <- rf_lag[(r_win + p0 + 1):(r_win + p0 + p)]

u_ha    <- numeric(3)

# Overall
u_ha[1] <- asset_alloc(actual_sim, rf_p, fc_ha_sim, fc_vol, gamma)$util

# Expansion
u_ha[2] <- asset_alloc(actual_sim[idx_exp], rf_p[idx_exp], fc_ha_sim[idx_exp], fc_vol[idx_exp], gamma)$util

# Recession
u_ha[3] <- asset_alloc(actual_sim[idx_rec], rf_p[idx_rec], fc_ha_sim[idx_rec], fc_vol[idx_rec], gamma)$util

# col order: KS | SIC | POOL-AVG | POOL-DMSFE | DI | SOP
fc_sim_a <- cbind(fc_ks_sim, fc_sic_sim, fc_avg_sim, fc_dmsfe_sim, fc_di_sim, fc_sop_sim)
fc_sim_b <- cbind(fc_ks_sim_ct, fc_sic_sim_ct, fc_avg_sim_ct, fc_dmsfe_sim_ct, fc_di_sim_ct, fc_sop_sim_ct)

u_a <- matrix(0, nm, 3)
u_b <- matrix(0, nm, 3)

for (i in seq_len(nm)) {
  
  # Overall
  u_a[i, 1] <- asset_alloc(actual_sim, rf_p, fc_sim_a[, i], fc_vol, gamma)$util
  u_b[i, 1] <- asset_alloc(actual_sim, rf_p, fc_sim_b[, i], fc_vol, gamma)$util
  
  # Expansion
  u_a[i, 2] <- asset_alloc(actual_sim[idx_exp], rf_p[idx_exp], fc_sim_a[idx_exp, i], fc_vol[idx_exp], gamma)$util
  u_b[i, 2] <- asset_alloc(actual_sim[idx_exp], rf_p[idx_exp], fc_sim_b[idx_exp, i], fc_vol[idx_exp], gamma)$util
  
  # Recession
  u_a[i, 3] <- asset_alloc(actual_sim[idx_rec], rf_p[idx_rec], fc_sim_a[idx_rec, i], fc_vol[idx_rec], gamma)$util
  u_b[i, 3] <- asset_alloc(actual_sim[idx_rec], rf_p[idx_rec], fc_sim_b[idx_rec, i], fc_vol[idx_rec], gamma)$util
  
}

u_ha_mat <- matrix(u_ha, nrow = nm, ncol = 3, byrow = TRUE)
delta_a  <- 1200 * (u_a - u_ha_mat)
delta_b  <- 1200 * (u_b - u_ha_mat)

# Output Table 6.2
methods <- c("Kitchen sink", "SIC", "POOL-AVG", "POOL-DMSFE", "Diffusion index", "Sum-of-parts")

fmt_row <- function(r2, pv, de) {
  
  sprintf("%6.2f [%.2f] %6.2f", r2, pv, de)
  
}

hdr <- sprintf(
  "%-17s  %-24s  %-24s  %-24s",
  "Method",
  "       Overall       ",
  "      Expansion      ",
  "      Recession      "
)

sub <- sprintf(
  "%-17s  %-24s  %-24s  %-24s",
  "",
  "R2OS[p]    Delta(%)  ",
  "R2OS[p]    Delta(%)  ",
  "R2OS[p]    Delta(%)  "
)

cat("\nTable 6.2 Monthly Multiple Economic Variables, 1957:01-2010:12\n")
cat(hdr, "\n")
cat(sub, "\n")

cat("\nPanel A: Unrestricted forecasts\n")
for (j in seq_len(nm)) {
  
  cat(sprintf("%-17s  %s  %s  %s\n",
    methods[j],
    fmt_row(r2_a[j, 1], r2_a[j, 2], delta_a[j, 1]),
    fmt_row(r2_a[j, 3], r2_a[j, 4], delta_a[j, 2]),
    fmt_row(r2_a[j, 5], r2_a[j, 6], delta_a[j, 3])
  ))
  
}

cat("\nPanel B: Forecasts with non-negativity restrictions\n")
for (j in seq_len(nm)) {
  
  cat(sprintf("%-16s  %s  %s  %s\n",
    methods[j],
    fmt_row(r2_b[j, 1], r2_b[j, 2], delta_b[j, 1]),
    fmt_row(r2_b[j, 3], r2_b[j, 4], delta_b[j, 2]),
    fmt_row(r2_b[j, 5], r2_b[j, 6], delta_b[j, 3])
  ))
  
}

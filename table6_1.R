#!/usr/bin/env Rscript
# ======================================================== #
#
#                 Replication of Table 6.1
#
#                 Gabriel E. Cabrera-Guzmán
#                The University of Manchester
#
#                        Spring, 2026
#
#                https://gcabrerag.rbind.io
#
# ------------------------------ #
# email: gabriel.cabreraguzman@postgrad.manchester.ac.uk
# ======================================================== #

# Run refs/load.R first to generate data/rz2013_data.rds.
#
# Panel A: Unrestricted predictive regression forecasts
# Panel B: Campbell & Thompson (2008) sign restrictions
# Columns: Overall | Expansion | Recession
#          R2OS(%) [p-val]  Delta(annual %)
#
# R2OS computed with LOG equity premium (paper p.353).
# Utility gain (Delta) computed with SIMPLE equity premium (footnote 28).

# Load packages
library("tidyverse")

# Load auxiliary functions
source("R/asset-alloc.R")
source("R/clark-west.R")

ds       <- readRDS("data/rz2013_data.rds")
data_log <- ds$data_log # log eq premium + 14 predictors
data_sim <- ds$data_sim # simple eq premium + 14 predictors
rf_lag   <- ds$rf_lag
rec      <- ds$rec

# Setup
y_log <- data_log$eq_prem
y_sim <- data_sim$eq_prem
x  <- select(data_log, -eq_prem) # same 14 predictors for both
n  <- nrow(data_log)
nc <- ncol(x)

# Full-sample beta signs (log premium, for CT sign restriction)
bf <- vapply(seq_len(nc), function(j) {
  
  ols(y_log[2:n], cbind(x[1:(n - 1), j], 1))[1, 1]
  
}, numeric(1))

bf[5] <- 1  # SVAR sign override

r_win <- (1946 - 1926) * 12 + 1 # 241: initial estimation window
p0 <- (1956 - 1946) * 12 # 120: pre-evaluation burn-in
p  <- n - r_win - p0 # evaluation period length
np <- p0 + p

# Forecast loop (log and simple premium run together)
fc_ha_log <- numeric(np)
fc_ec_log <- matrix(0, np, nc)
fc_ct_log <- matrix(0, np, nc)
fc_ha_sim <- numeric(np)
fc_ec_sim <- matrix(0, np, nc)
fc_ct_sim <- matrix(0, np, nc)

for (i in seq_len(np)) {
  
  ri <- r_win + i - 1
  fc_ha_log[i] <- mean(y_log[1:ri])
  fc_ha_sim[i] <- mean(y_sim[1:ri])
  xt <- as.matrix(x[1:(ri - 1),])

  for (j in seq_len(nc)) {
    
    xj <- cbind(xt[, j], 1)
    xp <- c(x[ri, j], 1)

    # --- Log premium ---
    b_l <- ols(y_log[2:ri], xj)[, 1]
    fc_ec_log[i, j] <- xp %*% b_l
    same_l <- (bf[j] > 0 && b_l[1] > 0) || (bf[j] < 0 && b_l[1] < 0)
    fc_ct_log[i, j] <- if (same_l) fc_ec_log[i, j] else fc_ha_log[i]
    
    if (fc_ct_log[i, j] < 0) fc_ct_log[i, j] <- 0

    # --- Simple premium ---
    b_s <- ols(y_sim[2:ri], xj)[, 1]
    fc_ec_sim[i, j] <- xp %*% b_s
    same_s <- (bf[j] > 0 && b_s[1] > 0) || (bf[j] < 0 && b_s[1] < 0)
    fc_ct_sim[i, j] <- if (same_s) fc_ec_sim[i, j] else fc_ha_sim[i]
    
    if (fc_ct_sim[i, j] < 0) fc_ct_sim[i, j] <- 0
    
  }
  
}

# --- Trim to evaluation period (1957:01 – 2010:12) ---
idx_ev <- seq(p0 + 1, np)

actual_log <- y_log[(r_win + p0 + 1):n]
actual_sim <- y_sim[(r_win + p0 + 1):n]

fc_ha_log <- fc_ha_log[idx_ev]
fc_ec_log <- fc_ec_log[idx_ev,]
fc_ct_log <- fc_ct_log[idx_ev,]
fc_ha_sim <- fc_ha_sim[idx_ev]
fc_ec_sim <- fc_ec_sim[idx_ev,]
fc_ct_sim <- fc_ct_sim[idx_ev,]

idx_rec <- which(rec$recession == 1)
idx_exp <- which(rec$recession == 0)

# --- R2OS statistics (log premium) ---
e_ha <- actual_log - fc_ha_log
e_ec <- matrix(actual_log, nrow = length(actual_log), ncol = nc) - fc_ec_log
e_ct <- matrix(actual_log, nrow = length(actual_log), ncol = nc) - fc_ct_log

r2_ec <- matrix(0, nc, 6)
r2_ct <- matrix(0, nc, 6)

for (j in seq_len(nc)) {
  
  # Overall  
  r2_ec[j, 1:2] <- r2os_stat(e_ha, e_ec[, j], fc_ha_log, fc_ec_log[, j])
  r2_ct[j, 1:2] <- r2os_stat(e_ha, e_ct[, j], fc_ha_log, fc_ct_log[, j])
  
  # Expansion
  r2_ec[j, 3:4] <- r2os_stat(e_ha[idx_exp], e_ec[idx_exp, j], fc_ha_log[idx_exp], fc_ec_log[idx_exp, j])
  r2_ct[j, 3:4] <- r2os_stat(e_ha[idx_exp], e_ct[idx_exp, j], fc_ha_log[idx_exp], fc_ct_log[idx_exp, j])
  
  # Recession
  r2_ec[j, 5:6] <- r2os_stat(e_ha[idx_rec], e_ec[idx_rec, j], fc_ha_log[idx_rec], fc_ec_log[idx_rec, j])
  r2_ct[j, 5:6] <- r2os_stat(e_ha[idx_rec], e_ct[idx_rec, j], fc_ha_log[idx_rec], fc_ct_log[idx_rec, j])
  
}

# --- Asset allocation / utility gain (simple premium) ---
gamma   <- 5
win_vol <- 12 * 5
fc_vol  <- numeric(p)

for (t in seq_len(p)) {
  
  idx_v <- (r_win + p0 + t - win_vol):(r_win + p0 + t - 1)
  yv <- y_sim[idx_v]
  fc_vol[t] <- mean(yv ^ 2) - mean(yv) ^ 2
  
}

rf_p <- rf_lag[(r_win + p0 + 1):(r_win + p0 + p)]

u_ha <- numeric(3)
u_ha[1] <- asset_alloc(actual_sim, rf_p, fc_ha_sim, fc_vol, gamma)$util

u_ha[2] <- asset_alloc(actual_sim[idx_exp], rf_p[idx_exp], fc_ha_sim[idx_exp], fc_vol[idx_exp], gamma)$util
u_ha[3] <- asset_alloc(actual_sim[idx_rec], rf_p[idx_rec], fc_ha_sim[idx_rec], fc_vol[idx_rec], gamma)$util

u_ec <- matrix(0, nc, 3)
u_ct <- matrix(0, nc, 3)

for (i in seq_len(nc)) {
  
  # --- Overall ---
  u_ec[i, 1] <- asset_alloc(actual_sim, rf_p, fc_ec_sim[, i], fc_vol, gamma)$util
  u_ct[i, 1] <- asset_alloc(actual_sim, rf_p, fc_ct_sim[, i], fc_vol, gamma)$util
  
  # --- Expansion ---
  u_ec[i, 2] <- asset_alloc(actual_sim[idx_exp], rf_p[idx_exp], fc_ec_sim[idx_exp, i], fc_vol[idx_exp], gamma)$util
  u_ct[i, 2] <- asset_alloc(actual_sim[idx_exp], rf_p[idx_exp], fc_ct_sim[idx_exp, i], fc_vol[idx_exp], gamma)$util
  
  # --- Recession ---
  u_ec[i, 3] <- asset_alloc(actual_sim[idx_rec], rf_p[idx_rec], fc_ec_sim[idx_rec, i], fc_vol[idx_rec], gamma)$util
  u_ct[i, 3] <- asset_alloc(actual_sim[idx_rec], rf_p[idx_rec], fc_ct_sim[idx_rec, i], fc_vol[idx_rec], gamma)$util
  
}

u_ha_mat <- matrix(u_ha, nrow = nc, ncol = 3, byrow = TRUE)
delta_ec <- 1200 * (u_ec - u_ha_mat)
delta_ct <- 1200 * (u_ct - u_ha_mat)

# Output Table 6.1 
predictors <- c(
  "log(DP)", "log(DY)", "log(EP)", "log(DE)", "SVAR", "BM", "NTIS", "TBL", "LTY", "LTR", "TMS", "DFY", "DFR", "INFL"
)

fmt_row <- function(r2, pv, de) {
  
  sprintf("%6.2f [%.2f] %6.2f", r2, pv, de)
  
}

hdr <- sprintf(
  "%-10s  %-24s  %-24s  %-24s",
  "Variable",
  "       Overall       ",
  "      Expansion      ",
  "      Recession      "
)

sub <- sprintf(
  "%-10s  %-24s  %-24s  %-24s",
  "",
  "R2OS[p]    Delta(%)  ",
  "R2OS[p]    Delta(%)  ",
  "R2OS[p]    Delta(%)  "
)

cat("\nTable 6.1 Monthly Individual Economic Variables, 1957:01-2010:12\n")
cat(hdr, "\n")
cat(sub, "\n")

cat("\nPanel A: Unrestricted predictive regression forecasts\n")
for (j in seq_len(nc)) {
  
  cat(sprintf("%-10s  %s  %s  %s\n",
    predictors[j],
    fmt_row(r2_ec[j, 1], r2_ec[j, 2], delta_ec[j, 1]),
    fmt_row(r2_ec[j, 3], r2_ec[j, 4], delta_ec[j, 2]),
    fmt_row(r2_ec[j, 5], r2_ec[j, 6], delta_ec[j, 3])
  ))
  
}

cat("\nPanel B: Campbell & Thompson (2008) restrictions\n")
for (j in seq_len(nc)) {
  
  cat(sprintf("%-10s  %s  %s  %s\n",
    predictors[j],
    fmt_row(r2_ct[j, 1], r2_ct[j, 2], delta_ct[j, 1]),
    fmt_row(r2_ct[j, 3], r2_ct[j, 4], delta_ct[j, 2]),
    fmt_row(r2_ct[j, 5], r2_ct[j, 6], delta_ct[j, 3])
  ))
  
}

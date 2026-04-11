#'---
#' title: Figure 6.1 – Cumulative Differences in Squared Forecast Errors
#'         Monthly U.S. Equity Premium, Individual Economic Variables, 1957:01–2010:12
#'
#' Rapach & Zhou (2013), Chapter 6, Handbook of Economic Forecasting Vol. 2A
#'
#' Black lines: unrestricted predictive regression forecast (Panel A of Table 6.1)
#' Gray lines:  Campbell & Thompson (2008) sign-restricted forecast (Panel B)
#' Shaded bars: NBER-dated recessions
#'
#' CDSFE_{j,tau} = sum_{s=1}^{tau} (e_HA,s)^2 - sum_{s=1}^{tau} (e_j,s)^2   [eq. 47]
#' Positive = predictive regression beats historical average in MSFE up to tau.
#'
#' Run data_prep.R first to generate data/rz2013_data.rds.
#'---

# Read packages
library("tidyverse")

# Read auxiliary functions
source("R/asset-alloc.R")   
source("R/clark-west.R")    

ds       <- readRDS("data/rz2013_data.rds")
data_log <- ds$data_log   # log eq premium + 14 predictors
rf_lag   <- ds$rf_lag
rec      <- ds$rec        # recession indicator, length = evaluation period

# Setup (identical to table6_1.R)
y_log <- data_log$eq_prem
x     <- select(data_log, -eq_prem)
n     <- nrow(data_log)
nc    <- ncol(x)

# Full-sample beta signs for CT restrictions
bf <- vapply(seq_len(nc), function(j) {
  
  ols(y_log[2:n], cbind(x[1:(n - 1), j], 1))[1, 1]
  
}, numeric(1))

bf[5] <- 1  # SVAR sign override

r_win <- (1946 - 1926) * 12 + 1 # 241
p0 <- (1956 - 1946) * 12 # 120
p  <- n - r_win - p0 # evaluation period length (648)
np <- p0 + p

# ── Forecast loop ───────────────────────────────────────────────────────────────
fc_ha_log <- numeric(np)
fc_ec_log <- matrix(0, np, nc)
fc_ct_log <- matrix(0, np, nc)

for (i in seq_len(np)) {
  
  ri <- r_win + i - 1
  
  fc_ha_log[i] <- mean(y_log[1:ri])
  xt <- as.matrix(x[1:(ri - 1), ])

  for (j in seq_len(nc)) {
    
    xj  <- cbind(xt[, j], 1)
    xp  <- c(as.numeric(x[ri, j]), 1)
    b_l <- ols(y_log[2:ri], xj)[, 1]
    fc_ec_log[i, j] <- xp %*% b_l
    same_l <- (bf[j] > 0 && b_l[1] > 0) || (bf[j] < 0 && b_l[1] < 0)
    fc_ct_log[i, j] <- if (same_l) fc_ec_log[i, j] else fc_ha_log[i]
    
    if (fc_ct_log[i, j] < 0) fc_ct_log[i, j] <- 0
    
  }
  
}

# Trim to evaluation period (1957:01 – 2010:12) 
idx_ev     <- seq(p0 + 1, np)
actual_log <- y_log[(r_win + p0 + 1):n]
fc_ha_log  <- fc_ha_log[idx_ev]
fc_ec_log  <- fc_ec_log[idx_ev,]
fc_ct_log  <- fc_ct_log[idx_ev,]

# Cumulative differences in squared forecast errors (CDSFE) 
e_ha <- actual_log - fc_ha_log
e_ec <- matrix(actual_log, nrow = p, ncol = nc) - fc_ec_log
e_ct <- matrix(actual_log, nrow = p, ncol = nc) - fc_ct_log

# cdsfe[t, j] = cumulative (HA SFE − predictor j SFE) through period t
cdsfe_ec <- matrix(0, p, nc)
cdsfe_ct <- matrix(0, p, nc)

for (j in seq_len(nc)) {
  
  cdsfe_ec[, j] <- cumsum(e_ha ^ 2) - cumsum(e_ec[, j] ^ 2)
  cdsfe_ct[, j] <- cumsum(e_ha ^ 2) - cumsum(e_ct[, j] ^ 2)
  
}

# Dates and recession shading 
dates <- seq(as.Date("1957-01-01"), by = "month", length.out = p)

# Recession episodes: consecutive runs of rec$recession == 1
rec_flag <- rec$recession
rec_starts <- dates[which(diff(c(0L, rec_flag)) == 1L)]
rec_ends   <- dates[which(diff(c(rec_flag, 0L)) == -1L)]

# Panel labels and column-major factor levels for facet_wrap(dir = "v")
pred_names   <- c("log(DP)", "log(DY)", "log(EP)", "log(DE)", "SVAR", "BM", 
                   "NTIS", "TBL", "LTY", "LTR", "TMS", "DFY", "DFR", "INFL")
pred_letters <- paste0("(", letters[1:14], ")")
pred_labels  <- paste(pred_letters, pred_names)

col_major_idx <- seq(1, 14, 1)
fct_levels    <- pred_labels[col_major_idx]

# Long-format data frame
long_df <- bind_rows(
  tibble(
    date      = rep(dates, nc),
    cdsfe     = as.vector(cdsfe_ec),
    predictor = rep(pred_labels, each = p),
    type      = "Unrestricted"
  ),
  tibble(
    date      = rep(dates, nc),
    cdsfe     = as.vector(cdsfe_ct),
    predictor = rep(pred_labels, each = p),
    type      = "CT-restricted"
  )
) |>
  mutate(
    predictor = factor(predictor, levels = fct_levels),
    type      = factor(type, levels = c("Unrestricted", "CT-restricted"))
  )

# Recession rectangles (ymin/ymax = -Inf/Inf spans whatever the panel y-range is)
rec_df <- tibble(xmin = rec_starts, xmax = rec_ends)

# Plot
ggplot(long_df, aes(x = date, y = cdsfe, colour = type)) +
  geom_rect(
    data        = rec_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill        = "grey85",
    alpha       = 0.6
  ) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_line(linewidth = 0.45) +
  scale_colour_manual(
    values = c("Unrestricted" = "black", "CT-restricted" = "grey45")
  ) +
  scale_y_continuous(
    limits = c(-0.0275, 0.0275), 
    breaks = seq(-0.02, 0.02, 0.01)
  ) + 
  scale_x_date(
    breaks      = as.Date(paste0(seq(1965, 2005, by = 10), "-01-01")),
    date_labels = "%Y"
  ) +
  facet_wrap(~ predictor, ncol = 4, scales = "free_y") +
  labs(
    x       = NULL,
    y       = NULL,
    colour  = NULL,
    title   = "Figure 6.1  Cumulative Differences in Squared Forecast Errors, 1957:01-2010:12"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position  = "bottom",
    strip.background = element_blank(),
    strip.text       = element_text(size = 7),
    axis.text        = element_text(size = 6),
    panel.grid       = element_blank(),
    plot.title       = element_text(size = 8), 
    plot.caption = element_text(hjust = 0)
  )

ggsave("figures/figure6_1.pdf", width = 8, height = 5)
cat("Saved: figure6_1.pdf\n")

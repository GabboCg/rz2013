#'---
#' title: Data preparation – refs/ table scripts
#' Rapach & Zhou (2013), Handbook of Economic Forecasting Vol. 2A
#'
#' Loads the Goyal-Welch Excel file, renames all long column names to
#' short ones, builds every dataset needed by Table6_1.R and Table6_2.R,
#' and saves them to data/refs_datasets.rds.
#'
#' Outputs (list):
#'   data_log  – log equity premium + 14 predictors  (for R2OS)
#'   data_sim  – simple equity premium + 14 predictors (for Delta)
#'   sink_data – 12-predictor kitchen-sink dataset
#'   sop_data  – sum-of-the-parts variables
#'   rf_lag    – lagged risk-free rate
#'   rec       – NBER recession indicator
#'---

# Read packages
library("tidyverse")
library("janitor")
library("readxl")

# Read dataset
raw <- readxl::read_xls("data-raw/handbook_data.xls", sheet = 1)

# Rename ALL long Goyal-Welch column names to short names 
raw_short <- raw |>
  clean_names() |>
  filter(date_yyyymm >= "192611" & date_yyyymm <= "201012") %>% 
  (function(x) as.data.frame(apply(x, 2, as.numeric))) |>
  rename(
    date     = date_yyyymm,
    mkt_ret  = crsp_s_p_500_value_weighted_return_with_dividends,
    mkt_ex   = crsp_s_p_500_value_weighted_return_excluding_dividends,
    d12      = x12_month_moving_sum_of_s_p_500_dividends,
    sp500    = s_p_500_index,
    e12      = x12_month_moving_sum_of_s_p_500_earnings,
    nber_rec = nber_recession_dummies,
    nber_pk  = nber_recession_dummies_with_peak_included,
    SVAR     = monthly_sum_of_squared_daily_returns_on_s_p_500_index,
    BM       = djia_book_to_market_value_ratio,
    NTIS     = net_equity_expansion,
    TBL      = x3_month_treasury_bill_yield_secondary_market,
    LTY      = long_term_government_bond_yield,
    LTR      = long_term_government_bond_return,
    AAA      = moodys_aaa_rated_corporate_bond_yield,
    BAA      = moodys_baa_rated_corporate_bond_yield,
    CORPR    = long_term_corporate_bond_return,
    INFL     = cpi_all_urban_consumers_inflation_rate,
    rf       = risk_free_rate
  )

# Log equity premium dataset (R2OS uses log premium, paper p.353)
data_log <- raw_short |>
  mutate(
    eq_prem  = log(1 + mkt_ret) - log(1 + dplyr::lag(rf)),
    DP       = log(d12) - log(sp500),
    EP       = log(e12) - log(sp500),
    TMS      = LTY - TBL,
    DFY      = BAA - AAA,
    DFR      = CORPR - LTR,
    INFL_lag = dplyr::lag(INFL),
    DY       = log(d12) - log(dplyr::lag(sp500)),
    DE       = log(d12) - log(e12)
  ) |>
  select(eq_prem, DP, DY, EP, DE, SVAR, BM, NTIS, TBL, LTY, LTR, TMS, DFY, DFR, INFL_lag) |>
  na.omit()

# Simple equity premium dataset (asset allocation uses simple, fn.28) 
data_sim <- raw_short |>
  mutate(
    eq_prem  = mkt_ret - dplyr::lag(rf),
    DP       = log(d12) - log(sp500),
    EP       = log(e12) - log(sp500),
    TMS      = LTY - TBL,
    DFY      = BAA - AAA,
    DFR      = CORPR - LTR,
    INFL_lag = dplyr::lag(INFL),
    DY       = log(d12) - log(dplyr::lag(sp500)),
    DE       = log(d12) - log(e12)
  ) |>
  select(eq_prem, DP, DY, EP, DE, SVAR, BM, NTIS, TBL, LTY, LTR, TMS, DFY, DFR, INFL_lag) |>
  na.omit()

# Kitchen-sink dataset (12 predictors: all except DE and TMS) 
sink_data <- raw_short |>
  mutate(
    eq_prem  = mkt_ret - dplyr::lag(rf),
    DP       = log(d12) - log(sp500),
    EP       = log(e12) - log(sp500),
    DFY      = BAA - AAA,
    DFR      = CORPR - LTR,
    INFL_lag = dplyr::lag(INFL),
    DY       = log(d12) - log(dplyr::lag(sp500))
  ) |>
  select(eq_prem, DP, DY, EP, SVAR, BM, NTIS, TBL, LTY, LTR, DFY, DFR, INFL_lag) |>
  na.omit()

# Sum-of-the-parts (Ferreira & Santa-Clara 2011)  
# Two versions mirror the two MATLAB scripts:
#   _log : log-return version used for R2OS  (Forecasts_monthly_log.m)
#   _sim : simple-return version used for Delta (Forecasts_monthly.m)
#
# e_growth_log = log(E12) - log(E12_lag)         [log earnings growth]
# dp_sop_log   = log(1 + (1/12)*D12/SP500)       [log(1+D/P)]
# rf_sop_log   = log(1 + rf)                     [log risk-free rate]
#
# e_growth_sim = (E12 - E12_lag) / E12_lag       [fractional earnings growth]
# dp_sop_sim   = (1/12) * D12 / SP500            [monthly D/P]
# rf_sop_sim   = rf                              [simple risk-free rate]
sop_data <- raw_short |>
  mutate(
    e12_lag      = dplyr::lag(e12, 1),
    e_growth_log = log(e12) - log(e12_lag),
    dp_sop_log   = log(1 + (1 / 12) * d12 / sp500),
    rf_sop_log   = log(1 + rf),
    e_growth_sim = (e12 - e12_lag) / e12_lag,
    dp_sop_sim   = (1 / 12) * d12 / sp500,
    rf_sop_sim   = rf
  ) |>
  select(e_growth_log, dp_sop_log, rf_sop_log, e_growth_sim, dp_sop_sim, rf_sop_sim) |>
  na.omit()

# Lagged risk-free rate (aligned to evaluation period) 
# Drop the last row (2010:12) to get the rf known at forecast formation
rf_lag <- raw_short$rf[-nrow(raw_short)]

# NBER recession indicator (1957:01 onward, peak-inclusive) 
rec <- raw_short |>
  filter(!is.na(nber_pk), date >= 195701) |>
  select(date, recession = nber_pk)

# Save to RDS
saveRDS(
  list(
    data_log  = data_log,
    data_sim  = data_sim,
    sink_data = sink_data,
    sop_data  = sop_data,
    rf_lag    = rf_lag,
    rec       = rec
  ),
  file = "data/rz2013_data.rds"
)

message("Saved: data/refs_datasets.rds")
message("  data_log:  ", nrow(data_log),  " x ", ncol(data_log))
message("  data_sim:  ", nrow(data_sim),  " x ", ncol(data_sim))
message("  sink_data: ", nrow(sink_data), " x ", ncol(sink_data))
message("  sop_data:  ", nrow(sop_data),  " x ", ncol(sop_data))
message("  rf_lag:    ", length(rf_lag),  " obs")
message("  rec:       ", nrow(rec),       " rows (", sum(rec$recession), " recession months)")

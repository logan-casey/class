
# ================================================================
# DOBRIDGE (2016): Replication of Tables 4, 5, 6
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(fixest)
  library(lubridate)
})

TAU <- 0.35

# ---- 1. LOAD & DEDUPLICATE -------------------------------------
comp_raw <- read_csv("compustat.csv", show_col_types = FALSE)
names(comp_raw) <- tolower(names(comp_raw))

comp_raw <- comp_raw %>%
  mutate(gvkey = as.integer(gvkey),
         fyear  = as.integer(fyear)) %>%
  mutate(
    sort_datafmt = case_when(
      tolower(coalesce(datafmt, "")) == "indl" ~ 0L,
      TRUE ~ 1L
    ),
    sort_consol = case_when(
      tolower(coalesce(consol, "")) == "c" ~ 0L,
      TRUE ~ 1L
    )
  ) %>%
  arrange(gvkey, fyear, sort_datafmt, sort_consol, desc(coalesce(at, 0))) %>%
  group_by(gvkey, fyear) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-sort_datafmt, -sort_consol) %>%
  arrange(gvkey, fyear)

cat(sprintf("Rows after dedup: %d | Unique gvkey-fyear: %d\n",
            nrow(comp_raw),
            n_distinct(paste(comp_raw$gvkey, comp_raw$fyear))))

# ---- 2. FF48 ---------------------------------------------------
ff48_from_sic <- function(sic) {
  s <- as.integer(sic)
  case_when(
    (s>=100&s<=299)|(s>=700&s<=799)|(s>=910&s<=919)|s==2048 ~ 1L,
    (s>=2000&s<=2046)|(s>=2050&s<=2063)|(s>=2070&s<=2079)|
      (s>=2090&s<=2095)|(s>=2098&s<=2099) ~ 2L,
    (s>=2064&s<=2068)|s==2086|s==2087|s==2096|s==2097 ~ 3L,
    (s>=2080&s<=2085) ~ 4L,
    (s>=2100&s<=2199) ~ 5L,
    (s>=920&s<=999)|(s>=3650&s<=3652)|s==3732|s==3949|
      (s>=7800&s<=7833)|(s>=7840&s<=7841)|
      (s>=7900&s<=7997)|s==7999 ~ 6L,
    (s>=7812&s<=7819)|(s>=7820&s<=7823) ~ 7L,
    (s>=2700&s<=2749)|(s>=2770&s<=2771)|(s>=2780&s<=2799) ~ 8L,
    s==2047|(s>=2391&s<=2392)|(s>=2510&s<=2519)|(s>=2590&s<=2599)|
      (s>=2840&s<=2844)|(s>=3160&s<=3162)|(s>=3170&s<=3172)|
      (s>=3190&s<=3199)|s==3229|s==3260|(s>=3262&s<=3263)|s==3269|
      (s>=3230&s<=3231)|(s>=3630&s<=3639)|(s>=3750&s<=3751)|s==3800|
      (s>=3860&s<=3861)|(s>=3870&s<=3873)|s==3910|(s>=3914&s<=3915)|
      (s>=3960&s<=3962)|s==3991|s==3995 ~ 9L,
    (s>=2300&s<=2390)|(s>=3020&s<=3021)|(s>=3100&s<=3111)|
      (s>=3130&s<=3131)|(s>=3140&s<=3149)|(s>=3150&s<=3151)|
      (s>=3963&s<=3965) ~ 10L,
    (s>=8000&s<=8099) ~ 11L,
    s==3836|(s>=3840&s<=3851)|s==5047|s==5122 ~ 12L,
    (s>=2830&s<=2836) ~ 13L,
    (s>=2800&s<=2829)|(s>=2850&s<=2879)|(s>=2890&s<=2899) ~ 14L,
    s==3031|s==3041|(s>=3050&s<=3053)|(s>=3060&s<=3089) ~ 15L,
    (s>=2200&s<=2284)|(s>=2290&s<=2299)|
      (s>=2393&s<=2395)|(s>=2397&s<=2399) ~ 16L,
    (s>=800&s<=899)|(s>=2400&s<=2439)|(s>=2450&s<=2459)|
      (s>=2490&s<=2499)|(s>=2660&s<=2661)|(s>=2950&s<=2952)|
      (s>=3200&s<=3229)|(s>=3240&s<=3241)|(s>=3250&s<=3259)|
      s==3261|s==3264|(s>=3270&s<=3275)|(s>=3280&s<=3281)|
      (s>=3290&s<=3293)|(s>=3295&s<=3299)|(s>=3420&s<=3429)|
      (s>=3430&s<=3433)|(s>=3440&s<=3442)|s==3446|s==3448|s==3449|
      (s>=3460&s<=3469)|(s>=3490&s<=3499)|s==3996 ~ 17L,
    (s>=1500&s<=1511)|(s>=1520&s<=1542)|(s>=1600&s<=1799) ~ 18L,
    (s>=3300&s<=3316)|(s>=3320&s<=3325)|(s>=3330&s<=3341)|
      (s>=3350&s<=3357)|(s>=3360&s<=3379)|(s>=3390&s<=3399) ~ 19L,
    (s>=3410&s<=3412)|(s>=3443&s<=3444)|
      (s>=3460&s<=3462)|(s>=3490&s<=3499) ~ 20L,
    (s>=3510&s<=3536)|s==3538|(s>=3540&s<=3569)|
      (s>=3580&s<=3582)|(s>=3585&s<=3586)|
      s==3589|(s>=3590&s<=3599) ~ 21L,
    (s>=3600&s<=3621)|(s>=3623&s<=3629)|(s>=3640&s<=3647)|
      (s>=3648&s<=3649)|s==3660|(s>=3690&s<=3699)|
      (s>=3810&s<=3812) ~ 22L,
    (s>=3711&s<=3716)|(s>=3750&s<=3751)|(s>=5010&s<=5015)|
      (s>=5510&s<=5521)|(s>=5530&s<=5531)|(s>=5560&s<=5561)|
      (s>=5570&s<=5571)|(s>=5580&s<=5581)|
      (s>=5590&s<=5599)|(s>=7510&s<=7549) ~ 23L,
    (s>=3720&s<=3727)|s==3728|(s>=3730&s<=3731) ~ 24L,
    (s>=3730&s<=3731)|(s>=3740&s<=3743) ~ 25L,
    (s>=3760&s<=3769)|(s>=3780&s<=3799)|s==3812 ~ 26L,
    (s>=1040&s<=1049) ~ 27L,
    (s>=1000&s<=1039)|(s>=1050&s<=1090)|
      (s>=1094&s<=1094)|(s>=1096&s<=1099)|(s>=1400&s<=1499) ~ 28L,
    (s>=1200&s<=1299) ~ 29L,
    (s>=1300&s<=1399)|(s>=2900&s<=2912)|(s>=2990&s<=2999)|
      (s>=5170&s<=5172)|s==5990 ~ 30L,
    (s>=4900&s<=4942)|(s>=4950&s<=4991) ~ 31L,
    (s>=4800&s<=4899) ~ 32L,
    (s>=7020&s<=7021)|(s>=7040&s<=7041)|s==7080|
      (s>=7200&s<=7299)|s==7395|s==7500|(s>=7600&s<=7699) ~ 33L,
    (s>=7370&s<=7374)|(s>=7376&s<=7380)|(s>=7389&s<=7390)|
      s==7392|s==7394|s==7396|s==7397|s==7399|s==8742 ~ 34L,
    (s>=3570&s<=3579)|(s>=3680&s<=3689)|s==3695 ~ 35L,
    s==7371|s==7372|s==7373|s==7375 ~ 36L,
    s==3622|(s>=3661&s<=3679)|(s>=3825&s<=3827)|
      s==3829|(s>=3841&s<=3851) ~ 37L,
    (s>=3811&s<=3813)|(s>=3820&s<=3826)|s==3828|s==3990 ~ 38L,
    (s>=2440&s<=2449)|(s>=2520&s<=2549)|(s>=2600&s<=2659)|
      (s>=2670&s<=2699)|(s>=2750&s<=2769)|(s>=3220&s<=3221)|
      (s>=3990&s<=3999) ~ 39L,
    (s>=2650&s<=2661)|(s>=3410&s<=3412) ~ 40L,
    (s>=4000&s<=4013)|(s>=4040&s<=4049)|(s>=4100&s<=4199)|
      (s>=4200&s<=4299)|(s>=4400&s<=4499)|(s>=4600&s<=4699)|
      (s>=4700&s<=4799) ~ 41L,
    (s>=5000&s<=5099)|(s>=5100&s<=5199) ~ 42L,
    (s>=5200&s<=5299)|(s>=5600&s<=5699)|(s>=5900&s<=5999) ~ 43L,
    (s>=5800&s<=5829)|s==5890|(s>=7000&s<=7019)|
      (s>=7040&s<=7049)|s==7213 ~ 44L,
    (s>=6000&s<=6099)|(s>=6020&s<=6026)|
      (s>=6120&s<=6179)|(s>=6180&s<=6199) ~ 45L,
    (s>=6310&s<=6411) ~ 46L,
    (s>=6500&s<=6553)|(s>=6725&s<=6726)|s==6798 ~ 47L,
    TRUE ~ 48L
  )
}

# ---- 3. TAXABLE INCOME -----------------------------------------
comp <- comp_raw %>%
  mutate(
    pidom_imp  = coalesce(
      pidom,
      if_else(!is.na(pi), pi - coalesce(pifo, 0L), NA_real_)
    ),
    txdfed_use = coalesce(txdfed, txdi),
    xido_use   = coalesce(xido, 0),
    ti = pidom_imp - coalesce(txdfed_use, 0) / TAU +
      xido_use / (1 - TAU)
  )

# ---- 4. TAX RATE -----------------------------------------------
comp <- comp %>%
  mutate(
    txfed_use = coalesce(
      txfed,
      coalesce(txt, 0) - coalesce(txfo, 0)
      - coalesce(txs,  0)
      - coalesce(txo,  0)
    ),
    tax_rate = case_when(
      !is.na(ti) & ti > 0 &
        !is.na(txfed_use) & txfed_use > 0 ~ txfed_use / ti,
      TRUE ~ NA_real_
    ),
    tax_rate = if_else(tax_rate < 0.01 | tax_rate > 0.52, NA_real_, tax_rate)
  )

# ---- 5. MARGINAL TAX RATE (Graham & Mills 2008) ----------------
comp <- comp %>%
  mutate(
    tax_ratio_gm = case_when(
      !is.na(txfed_use) & !is.na(pidom_imp) & pidom_imp != 0
      ~ txfed_use / pidom_imp,
      !is.na(txt) & !is.na(pi) & pi != 0
      ~ txt / pi,
      TRUE ~ NA_real_
    ),
    is_multinational = !is.na(pi) & pi != 0 &
      abs(coalesce(pifo, 0) / pi) > 0.05,
    mtr = 0.331 -
      0.075 * as.integer(!is.na(tax_ratio_gm) & tax_ratio_gm < 0.1) -
      0.012 * as.integer(!is.na(tlcf) & tlcf > 0) -
      0.106 * as.integer(!is.na(pi)   & pi   < 0) +
      0.037 * as.integer(is_multinational)
  )

# ---- 6. OUTCOMES & CONTROLS ------------------------------------
comp <- comp %>%
  group_by(gvkey) %>%
  mutate(
    che_lag   = lag(che),
    dltt_lag  = lag(coalesce(dltt, 0)),
    dlc_lag   = lag(coalesce(dlc,  0)),
    ivst_lag  = lag(coalesce(ivst, 0)),
    ppent_lag = lag(ppent),
    emp_lag   = lag(emp)
  ) %>%
  ungroup() %>%
  mutate(
    capxv_use        = coalesce(capxv, capx, 0),
    investment       = capxv_use - coalesce(sppe, 0),
    delta_cash       = coalesce(che, 0) - coalesce(che_lag, 0),
    dltt_use         = coalesce(dltt, 0),
    dlc_use          = coalesce(dlc,  0),
    delta_debt_total = (dltt_use + dlc_use) - (dltt_lag + dlc_lag),
    delta_debt_lt    = dltt_use - dltt_lag,
    delta_debt_st    = dlc_use  - dlc_lag,
    payout           = coalesce(dvc, 0) + coalesce(prstkc, 0),
    delta_st_invest  = coalesce(ivst, 0) - ivst_lag,
    delta_lt_invest  = coalesce(ivch, 0),
    acquisitions     = coalesce(aqc,  0),
    other_uses       = acquisitions +
      coalesce(delta_st_invest, 0) +
      coalesce(delta_lt_invest, 0),
    delta_emp        = emp - emp_lag,
    at_use           = coalesce(at,     NA_real_),
    prccf_use        = coalesce(prcc_f, prcc_c, NA_real_),
    seq_use          = coalesce(seq,    NA_real_),
    tobins_q = (at_use + prccf_use * coalesce(csho, NA_real_) -
                  (seq_use + coalesce(txditc, 0) - coalesce(pstk, 0))) / at_use,
    roa          = coalesce(oibdp, 0) / at_use,
    cf_assets    = (coalesce(ib, 0) + coalesce(dp, 0)) / at_use,
    sales_assets = coalesce(sale, 0) / at_use,
    leverage     = (dlc_use + dltt_use) /
      (dlc_use + dltt_use + prccf_use * coalesce(csho, NA_real_)),
    ln_at        = log(at_use),
    cash_assets  = coalesce(che, 0) / at_use,
    ppe_assets   = coalesce(ppent, 0) / at_use,
    kz_index = case_when(
      !is.na(ppent_lag) & ppent_lag > 0 ~ (
        -1.002 * (coalesce(ib,0) + coalesce(dp,0)) / ppent_lag +
          0.283 * (at_use + prccf_use * coalesce(csho,NA_real_) -
                     coalesce(ceq,0) - coalesce(txdb,0)) / at_use +
          3.139 * (dltt_use + dlc_use) /
          (dltt_use + coalesce(seq_use,0)) -
          39.368 * (coalesce(dvc,0) + coalesce(dvp,0)) / ppent_lag -
          1.315 * coalesce(che,0) / ppent_lag
      ),
      TRUE ~ NA_real_
    ),
    sic_use = coalesce(as.integer(sich), as.integer(sic), NA_integer_),
    ff48    = ff48_from_sic(sic_use)
  )

# ---- 7. REFUND CALCULATION FUNCTIONS ---------------------------
adjust_within_window_2yr <- function(ti_vec) {
  ti_vec <- replace(ti_vec, is.na(ti_vec), 0)
  avail  <- pmax(ti_vec, 0)
  for (j in seq_along(ti_vec)) {
    if (ti_vec[j] >= 0) next
    loss_rem <- -ti_vec[j]
    if (j >= 3 && avail[j-2] > 0) {
      use        <- min(avail[j-2], loss_rem)
      avail[j-2] <- avail[j-2] - use
      loss_rem   <- loss_rem - use
    }
    if (j >= 2 && loss_rem > 0 && avail[j-1] > 0) {
      use        <- min(avail[j-1], loss_rem)
      avail[j-1] <- avail[j-1] - use
    }
  }
  avail
}

calc_refund_and_V <- function(ti_vec, tr_vec, policy_loss, cap_first = 1) {
  if (is.na(policy_loss) || policy_loss <= 0)
    return(list(refund = 0, V = NA_real_))
  ti_vec <- replace(ti_vec, is.na(ti_vec), 0)
  tr_vec <- replace(tr_vec, is.na(tr_vec), TAU)
  tr_vec <- pmin(pmax(tr_vec, 0), 0.52)
  avail    <- adjust_within_window_2yr(ti_vec)
  avail[1] <- avail[1] * cap_first
  V        <- sum(avail) - policy_loss
  refund   <- 0
  remaining <- policy_loss
  for (i in seq_along(avail)) {
    if (remaining < 1e-9) break
    use       <- min(avail[i], remaining)
    refund    <- refund + use * tr_vec[i]
    remaining <- remaining - use
  }
  list(refund = refund, V = V)
}

compute_usage <- function(ti_vec, policy_loss) {
  policy_loss <- if (is.na(policy_loss) || policy_loss < 0) 0 else policy_loss
  ti_vec      <- replace(ti_vec, is.na(ti_vec), 0)
  avail       <- pmax(adjust_within_window_2yr(ti_vec), 0)
  remaining   <- policy_loss
  usage       <- numeric(length(ti_vec))
  for (i in seq_along(avail)) {
    if (remaining < 1e-9) break
    use      <- min(avail[i], remaining)
    usage[i] <- use
    remaining <- remaining - use
  }
  usage
}

# ---- 8. DATASET CONSTRUCTION HELPERS ---------------------------
comp_thin <- comp %>%
  dplyr::select(
    gvkey, fyear, ti, tax_rate, mtr, at,
    tobins_q, roa, cf_assets, sales_assets, leverage, ln_at,
    cash_assets, ppe_assets, kz_index, is_multinational, ff48, sic_use,
    investment, delta_cash, delta_debt_total, delta_debt_lt, delta_debt_st,
    payout, other_uses, delta_emp, dvc, prstkc, oibdp
  )

stopifnot(
  "comp_thin still has duplicate gvkey/fyear — check dedup step" =
    !anyDuplicated(paste(comp_thin$gvkey, comp_thin$fyear)))

get_window_wide <- function(df, window_yrs) {
  df %>%
    filter(fyear %in% window_yrs) %>%
    dplyr::select(gvkey, fyear, ti, mtr) %>%
    group_by(gvkey, fyear) %>%
    slice(1) %>%
    ungroup() %>%
    pivot_wider(
      names_from  = fyear,
      values_from = c(ti, mtr),
      names_glue  = "{.value}_{fyear}"
    )
}

apply_refund_calc <- function(wide_df, window_yrs, cap_first = 1) {
  ti_cols  <- paste0("ti_",  window_yrs)
  mtr_cols <- paste0("mtr_", window_yrs)
  n        <- length(window_yrs)
  for (col in c(ti_cols, mtr_cols))
    if (!col %in% names(wide_df)) wide_df[[col]] <- NA_real_
  input_mat <- wide_df %>%
    dplyr::select(all_of(c(ti_cols, mtr_cols, "policy_loss"))) %>%
    as.data.frame()
  results <- pmap(input_mat, function(...) {
    vals        <- list(...)
    ti_vec      <- as.numeric(unlist(vals[seq_len(n)],     use.names = FALSE))
    mtr_vec     <- as.numeric(unlist(vals[seq_len(n) + n], use.names = FALSE))
    policy_loss <- as.numeric(vals[[2 * n + 1]])
    calc_refund_and_V(ti_vec, mtr_vec, policy_loss, cap_first)
  })
  wide_df %>%
    mutate(
      TaxRefund = map_dbl(results, "refund"),
      V         = map_dbl(results, "V")
    )
}

merge_outcomes <- function(calc_df, loss_yr, outcome_yr) {
  outcomes <- comp_thin %>%
    filter(fyear == outcome_yr) %>%
    dplyr::select(gvkey, investment, delta_cash, delta_debt_total,
                  delta_debt_lt, delta_debt_st, payout, other_uses,
                  delta_emp, at)
  controls <- comp_thin %>%
    filter(fyear == loss_yr) %>%
    dplyr::select(gvkey,
                  tobins_q_pre     = tobins_q,
                  roa_pre          = roa,
                  cf_assets_pre    = cf_assets,
                  sales_assets_pre = sales_assets,
                  leverage_pre     = leverage,
                  mtr_pre          = mtr,
                  ln_at_pre        = ln_at,
                  cash_assets_pre  = cash_assets,
                  ppe_assets_pre   = ppe_assets,
                  kz_pre           = kz_index,
                  dvc_pre          = dvc,
                  prstkc_pre       = prstkc,
                  oibdp_pre        = oibdp,
                  is_multi         = is_multinational,
                  ff48, sic_use)
  calc_df %>%
    dplyr::select(gvkey, TaxRefund, V, policy_loss) %>%
    mutate(losses_sq    = policy_loss^2,
           loss_year    = as.integer(loss_yr),
           outcome_year = as.integer(outcome_yr)) %>%
    inner_join(outcomes, by = "gvkey") %>%
    inner_join(controls, by = "gvkey") %>%
    filter(!is.na(at), at > 1)
}

# ---- 8A. 2002 POLICY PERIOD ------------------------------------
loss_2001 <- comp_thin %>%
  filter(fyear == 2001, !is.na(ti), ti < 0) %>%
  transmute(gvkey, policy_loss = -ti)

df_2001_calc <- loss_2001 %>%
  inner_join(get_window_wide(comp_thin, 1996:2000), by = "gvkey") %>%
  filter(if_all(paste0("ti_", 1996:2000), ~ !is.na(.))) %>%
  apply_refund_calc(1996:2000, cap_first = 1)

# Define loss_2002 before §3C adjustment
loss_2002 <- comp_thin %>%
  filter(fyear == 2002, !is.na(ti), ti < 0) %>%
  transmute(gvkey, policy_loss = -ti)

# §3C cross-year adjustment: dual-loss firms only
dual_loss_gvkeys <- intersect(loss_2001$gvkey, loss_2002$gvkey)
cat(sprintf("Firms with losses in both 2001 and 2002: %d\n",
            length(dual_loss_gvkeys)))

consumed_by_2001 <- df_2001_calc %>%
  filter(gvkey %in% dual_loss_gvkeys) %>%
  dplyr::select(gvkey, policy_loss,
                ti_1996, ti_1997, ti_1998, ti_1999, ti_2000) %>%
  as.data.frame() %>%
  pmap_dfr(function(gvkey, policy_loss,
                    ti_1996, ti_1997, ti_1998, ti_1999, ti_2000) {
    ti_vec <- c(ti_1996, ti_1997, ti_1998, ti_1999, ti_2000)
    usage  <- compute_usage(ti_vec, policy_loss)
    tibble(gvkey         = gvkey,
           consumed_1997 = usage[2],
           consumed_1998 = usage[3],
           consumed_1999 = usage[4],
           consumed_2000 = usage[5])
  })

cat(sprintf("Consumed_by_2001 rows: %d\n", nrow(consumed_by_2001)))

win_2002_adj <- get_window_wide(comp_thin, 1997:2001) %>%
  left_join(consumed_by_2001, by = "gvkey") %>%
  mutate(
    ti_1997 = pmax(coalesce(ti_1997, 0) - coalesce(consumed_1997, 0), 0),
    ti_1998 = pmax(coalesce(ti_1998, 0) - coalesce(consumed_1998, 0), 0),
    ti_1999 = pmax(coalesce(ti_1999, 0) - coalesce(consumed_1999, 0), 0),
    ti_2000 = pmax(coalesce(ti_2000, 0) - coalesce(consumed_2000, 0), 0),
    ti_2001 = pmax(coalesce(ti_2001, 0), 0)
  ) %>%
  dplyr::select(-starts_with("consumed_"))

df_2002_calc <- loss_2002 %>%
  inner_join(win_2002_adj, by = "gvkey") %>%
  filter(if_all(paste0("ti_", 1997:2001), ~ !is.na(.))) %>%
  apply_refund_calc(1997:2001, cap_first = 1)

df_2002_full <- bind_rows(
  merge_outcomes(df_2001_calc, 2001L, 2002L),
  merge_outcomes(df_2002_calc, 2002L, 2003L)
)

# ---- 8B. 2009 POLICY PERIOD ------------------------------------
make_2009_candidates <- function(loss_yr, win_yrs) {
  comp_thin %>%
    filter(fyear == loss_yr, !is.na(ti), ti < 0) %>%
    transmute(gvkey, policy_loss = -ti) %>%
    inner_join(get_window_wide(comp_thin, win_yrs), by = "gvkey") %>%
    filter(if_all(paste0("ti_", win_yrs), ~ !is.na(.))) %>%
    apply_refund_calc(win_yrs, cap_first = 0.5) %>%
    dplyr::select(gvkey, TaxRefund, V, policy_loss)
}

df_2009_full <- bind_rows(
  make_2009_candidates(2008, 2003:2007),
  make_2009_candidates(2009, 2004:2008)
) %>%
  group_by(gvkey) %>%
  slice_max(TaxRefund, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(losses_sq    = policy_loss^2,
         loss_year    = NA_integer_,
         outcome_year = 2010L) %>%
  inner_join(
    comp_thin %>%
      filter(fyear == 2010) %>%
      dplyr::select(gvkey, investment, delta_cash, delta_debt_total,
                    delta_debt_lt, delta_debt_st, payout, other_uses,
                    delta_emp, at),
    by = "gvkey"
  ) %>%
  inner_join(
    comp_thin %>%
      filter(fyear == 2009) %>%
      dplyr::select(gvkey,
                    tobins_q_pre     = tobins_q,
                    roa_pre          = roa,
                    cf_assets_pre    = cf_assets,
                    sales_assets_pre = sales_assets,
                    leverage_pre     = leverage,
                    mtr_pre          = mtr,
                    ln_at_pre        = ln_at,
                    cash_assets_pre  = cash_assets,
                    ppe_assets_pre   = ppe_assets,
                    kz_pre           = kz_index,
                    dvc_pre          = dvc,
                    prstkc_pre       = prstkc,
                    oibdp_pre        = oibdp,
                    is_multi         = is_multinational,
                    ff48, sic_use),
    by = "gvkey"
  ) %>%
  filter(!is.na(at), at > 1)

cat(sprintf("2009 period (pre-trim): %d firms\n", nrow(df_2009_full)))

# ---- 8C. CRSP MERGE --------------------------------------------
# Require active CRSP listing in the outcome year.
# Load crosswalk
xwalk <- read_csv("crsp_compustat_xwalk.csv", show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  mutate(gvkey = as.integer(gvkey)) %>%
  filter(linktype %in% c("LC","LU","LS"),
         linkprim %in% c("P","C")) %>%
  mutate(
    linkdt    = as.Date(linkdt),
    linkenddt = if_else(is.na(linkenddt) | linkenddt == as.Date("2099-12-31"),
                        Sys.Date(), as.Date(linkenddt))
  ) %>%
  dplyr::select(gvkey, lpermno, linkdt, linkenddt)

# Expand crosswalk to gvkey-year panel
crsp_gvkeys_by_year <- xwalk %>%
  mutate(
    start_yr = year(linkdt),
    end_yr   = case_when(
      is.na(linkenddt)       ~ year(Sys.Date()),
      year(linkenddt) > 2030 ~ year(Sys.Date()),
      TRUE                   ~ year(linkenddt)
    )
  ) %>%
  filter(!is.na(start_yr), !is.na(end_yr), end_yr >= start_yr) %>%
  rowwise() %>%
  mutate(fyear = list(seq(start_yr, end_yr))) %>%
  ungroup() %>%
  unnest(fyear) %>%
  group_by(gvkey, fyear) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(gvkey, fyear, lpermno)

cat(sprintf("xwalk gvkey-years: %d\n", nrow(crsp_gvkeys_by_year)))

# Load CRSP monthly stock file
crsp_msf_raw <- read_csv("crsp_msf.csv", show_col_types = FALSE)

crsp_msf <- crsp_msf_raw %>%
  mutate(
    permno = as.integer(PERMNO),
    fyear  = year(date),
    prc    = PRC
  ) %>%
  dplyr::select(permno, fyear, prc)

# Active = has at least one valid price in that year
crsp_active_by_year <- crsp_msf %>%
  filter(!is.na(prc), abs(prc) > 0) %>%
  distinct(permno, fyear)

# Merge active listings into crosswalk
crsp_gvkeys_active <- crsp_gvkeys_by_year %>%
  inner_join(crsp_active_by_year,
             by = c("lpermno" = "permno", "fyear"))

cat(sprintf("Total gvkey-years in xwalk: %d\n", nrow(crsp_gvkeys_by_year)))
cat(sprintf("Active gvkey-years (CRSP):  %d\n", nrow(crsp_gvkeys_active)))

# Filter: require active CRSP listing in outcome year only
filter_to_crsp <- function(df) {
  crsp_outcome <- crsp_gvkeys_active %>%
    rename(outcome_year = fyear) %>%
    dplyr::select(gvkey, outcome_year)
  df %>%
    inner_join(crsp_outcome, by = c("gvkey", "outcome_year"))
}

df_2002_full_crsp <- filter_to_crsp(df_2002_full)
df_2009_full_crsp <- filter_to_crsp(
  df_2009_full %>% mutate(outcome_year = 2010L)
)

cat(sprintf("After CRSP filter — 2002: %d (was %d) | 2009: %d (was %d)\n",
            nrow(df_2002_full_crsp), nrow(df_2002_full),
            nrow(df_2009_full_crsp), nrow(df_2009_full)))

# ---- 9. TRIM 1% TAILS & RKD VARIABLES -------------------------
trim_and_prep <- function(df) {
  qv <- quantile(df$V,          c(.01,.99), na.rm=TRUE)
  qi <- quantile(df$investment, c(.01,.99), na.rm=TRUE)
  df %>%
    filter(between(V,          qv[1], qv[2]),
           between(investment, qi[1], qi[2])) %>%
    mutate(
      D    = as.integer(V < 0),
      V2   = V^2,
      D_V  = D * V,
      D_V2 = D * V^2
    )
}

df_2002 <- trim_and_prep(df_2002_full_crsp)
df_2009 <- trim_and_prep(df_2009_full_crsp)

cat(sprintf("Final — 2002: %d | 2009: %d\n", nrow(df_2002), nrow(df_2009)))

df_2002 %>%
  group_by(loss_year) %>%
  summarise(n=n(), frac_D=mean(D), mean_refund=mean(TaxRefund),
            mean_V=mean(V), mean_profits=mean(V+policy_loss)) %>%
  print()

# ---- 10. REGRESSION SETUP --------------------------------------
ctrl <- paste(
  "tobins_q_pre + roa_pre + cf_assets_pre + sales_assets_pre",
  "leverage_pre + mtr_pre + ln_at_pre + policy_loss + losses_sq",
  sep = " + "
)

run_iv <- function(outcome, df, year_fe = TRUE) {
  # Winsorize non-investment outcomes (replace outliers, don't drop rows)
  if (!outcome %in% c("investment", "delta_emp") && outcome %in% names(df)) {
    q_out <- quantile(df[[outcome]], c(.01,.99), na.rm=TRUE)
    df <- df %>%
      mutate(across(all_of(outcome),
                    ~pmin(pmax(., q_out[1]), q_out[2])))
  }
  fe_str <- if (year_fe && n_distinct(df$outcome_year) > 1)
    "ff48 + outcome_year"
  else
    "ff48"
  fml <- as.formula(
    paste0(outcome, " ~ V + V2 + ", ctrl,
           " | ", fe_str, " | TaxRefund ~ D_V + D_V2")
  )
  tryCatch(
    feols(fml, data = df, cluster = ~ff48),
    error = function(e) {
      message("  [IV failed: '", outcome, "'] ", conditionMessage(e))
      NULL
    }
  )
}

get_iv_coef <- function(m) {
  if (is.null(m)) return(c(coef=NA, se=NA))
  cn <- names(coef(m))
  nm <- cn[cn %in% c("TaxRefund", "fit_TaxRefund")]
  if (length(nm) == 0) {
    message("Available coef names: ", paste(cn, collapse=", "))
    return(c(coef=NA, se=NA))
  }
  c(coef = coef(m)[[nm[1]]], se = se(m)[[nm[1]]])
}

extract_row <- function(m, label) {
  if (is.null(m))
    return(tibble(label=label, coef=NA_real_, se=NA_real_,
                  sig="", n=NA_integer_, rsq=NA_real_))
  cn  <- names(coef(m))
  nm  <- cn[cn %in% c("TaxRefund", "fit_TaxRefund")]
  if (length(nm) == 0) {
    message("  [no IV coef in: ", paste(cn, collapse=", "), "]")
    return(tibble(label=label, coef=NA_real_, se=NA_real_,
                  sig="[no coef]", n=m$nobs, rsq=NA_real_))
  }
  nm <- nm[1]
  p  <- tryCatch(pvalue(m)[[nm]], error=function(e) NA_real_)
  sg <- case_when(is.na(p)~"", p<.01~"***", p<.05~"**", p<.10~"*", TRUE~"")
  tibble(label=label, coef=coef(m)[[nm]], se=se(m)[[nm]],
         sig=sg, n=m$nobs, rsq=r2(m, type="r2"))
}

print_panel <- function(rows, title) {
  cat("\n", strrep("─", 72), "\n", sep="")
  cat(title, "\n")
  cat(strrep("─", 72), "\n")
  cat(sprintf("%-32s %8s %8s %4s %7s %6s\n",
              "Outcome / Subsample", "Coeff", "SE", "Sig", "N", "R²"))
  cat(strrep("─", 72), "\n")
  pwalk(rows, function(label, coef, se, sig, n, rsq) {
    if (is.na(coef))
      cat(sprintf("%-32s  [failed]\n", label))
    else
      cat(sprintf("%-32s %8.3f %8.3f %4s %7d %6.3f\n",
                  label, coef, se, sig, n, rsq))
  })
}

# ---- 11. TABLE 4 — FIRST STAGE ---------------------------------
fs_fml_2002 <- as.formula(
  paste0("TaxRefund ~ V + V2 + D_V + D_V2 + ", ctrl, " | ff48 + outcome_year"))
fs_fml_2009 <- as.formula(
  paste0("TaxRefund ~ V + V2 + D_V + D_V2 + ", ctrl, " | ff48"))

fs_2002 <- feols(fs_fml_2002, data = df_2002, cluster = ~ff48)
fs_2009 <- feols(fs_fml_2009, data = df_2009, cluster = ~ff48)

run_iv_2009 <- function(outcome, df) run_iv(outcome, df, year_fe = FALSE)

m_iv_2002 <- run_iv("investment", df_2002)
m_iv_2009 <- run_iv("investment", df_2009)

cat("\nFirst-stage F-statistics for excluded instruments (should be >> 10):\n")
cat("2002 period:", fitstat(m_iv_2002, "ivf")[[1]]$stat, "\n")
cat("2009 period:", fitstat(m_iv_2009, "ivf")[[1]]$stat, "\n")

cat("\n", strrep("═", 72), "\n", sep="")
cat("TABLE 4: FIRST-STAGE — AVERAGE FIRM TAX RATE\n")
cat(strrep("═", 72), "\n")
cat(sprintf("%-22s %9s %9s %4s %7s %6s\n",
            "Period","D_V coef","SE","Sig","N","R²"))
cat(strrep("─", 72), "\n")

for (nm in c("2002","2009")) {
  m  <- if (nm == "2002") fs_2002 else fs_2009
  cn <- names(coef(m))
  if (!"D_V" %in% cn) {
    cat(sprintf("%-22s  [D_V not found]\n", paste(nm,"Policy"))); next
  }
  p  <- tryCatch(pvalue(m)[["D_V"]], error = function(e) NA_real_)
  sg <- case_when(is.na(p)~"", p<.01~"***", p<.05~"**", p<.10~"*", TRUE~"")
  cat(sprintf("%-22s %9.3f %9.3f %4s %7d %6.3f\n",
              paste(nm, "Policy"),
              coef(m)["D_V"], se(m)["D_V"],
              sg, m$nobs, r2(m, type="r2")))
}
cat("Paper: 0.337*** (2002) | 0.309*** (2009)\n")

# ---- 12. TABLE 5 — 2002 PERIOD ALLOCATION ----------------------
t5_vars <- c(
  "Investment"   = "investment",
  "Δ Cash"       = "delta_cash",
  "Δ Total Debt" = "delta_debt_total",
  "Payout"       = "payout",
  "Other Uses"   = "other_uses",
  "Δ Employment" = "delta_emp"
)

t5_rows <- imap_dfr(
  lapply(t5_vars, run_iv, df = df_2002, year_fe = TRUE),
  extract_row
)

print_panel(t5_rows,
            "TABLE 5: TAX REFUND ALLOCATION — 2002 POLICY PERIOD (refund years 2002-2003)")
cat("Paper: Investment 0.403*** | all others n.s.\n")

# ---- 13. TABLE 6 — INVESTMENT BY SUBGROUP ----------------------
df_2002 <- df_2002 %>%
  mutate(
    payout_pre = coalesce(dvc_pre, 0) + coalesce(prstkc_pre, 0),
    dff_ratio  = if_else(!is.na(oibdp_pre) & oibdp_pre != 0,
                         payout_pre / abs(oibdp_pre), NA_real_)
  )

kz_q75   <- quantile(df_2002$kz_pre,       .75, na.rm = TRUE)
kz_q25   <- quantile(df_2002$kz_pre,       .25, na.rm = TRUE)
tobq_q75 <- quantile(df_2002$tobins_q_pre, .75, na.rm = TRUE)
tobq_q25 <- quantile(df_2002$tobins_q_pre, .25, na.rm = TRUE)

subsets_a <- list(
  "KZ: constrained (top Q)"    = df_2002 %>% filter(kz_pre >= kz_q75),
  "KZ: unconstrained (bot Q)"  = df_2002 %>% filter(kz_pre <= kz_q25),
  "Payout: zero (constrained)" = df_2002 %>% filter(payout_pre == 0),
  "Payout: >0 (unconstrained)" = df_2002 %>% filter(payout_pre > 0),
  "DFF: ≤0 (constrained)"      = df_2002 %>%
    filter(!is.na(dff_ratio) & dff_ratio <= 0),
  "DFF: >0 (unconstrained)"    = df_2002 %>%
    filter(!is.na(dff_ratio) & dff_ratio >  0)
)

subsets_b <- list(
  "High Tobin's Q (top Q)"     = df_2002 %>% filter(tobins_q_pre >= tobq_q75),
  "Low Tobin's Q (bot Q)"      = df_2002 %>% filter(tobins_q_pre <= tobq_q25),
  "Refund year: 2002"          = df_2002 %>% filter(outcome_year == 2002),
  "Refund year: 2003"          = df_2002 %>% filter(outcome_year == 2003),
  "Multinational"              = df_2002 %>% filter(is_multi == TRUE),
  "Domestic"                   = df_2002 %>% filter(is_multi == FALSE)
)

rows_a <- imap_dfr(lapply(subsets_a, run_iv, outcome = "investment",
                          year_fe = TRUE), extract_row)
rows_b <- imap_dfr(lapply(subsets_b, run_iv, outcome = "investment",
                          year_fe = TRUE), extract_row)

cat("\n", strrep("═", 72), "\n", sep="")
cat("TABLE 6: INVESTMENT BY SUBGROUP — 2002 POLICY PERIOD\n")
print_panel(rows_a, "PANEL A — Financial Constraints")
cat("Paper: constrained $1.00-1.09* to **  |  unconstrained n.s.\n")
print_panel(rows_b, "PANEL B — Investment Opportunities")
cat("Paper: high-Q $1.02* | 2003 $0.75*** | multinational $0.68* | rest n.s.\n")
cat("\n", strrep("═", 72), "\n", "Done.\n", sep="")

# ---- 14. FIGURE 2 & FIGURE 3 -----------------------------------
library(ggplot2)

# Fixed residualize: pure zero-centered residuals, no grand mean added back
residualize <- function(var, df, year_fe = TRUE) {
  fe_str <- if (year_fe && n_distinct(df$outcome_year) > 1)
    "ff48 + outcome_year" else "ff48"
  fml    <- as.formula(paste0(var, " ~ ", ctrl, " | ", fe_str))
  m      <- feols(fml, data = df)
  rhs_vars <- c("tobins_q_pre","roa_pre","cf_assets_pre","sales_assets_pre",
                "leverage_pre","mtr_pre","ln_at_pre","policy_loss","losses_sq")
  df_out   <- df[!is.na(df[[var]]) & complete.cases(df[, rhs_vars]), ]
  df_out$resid_y <- resid(m)   # pure residuals, mean zero
  df_out[, c("V", "resid_y")]
}

# Require at least 10 obs per bin to suppress exploding CIs
make_fig2_bins <- function(resid_df, bin_width = 1, clip = 50, min_n = 10) {
  resid_df %>%
    filter(!is.na(V), !is.na(resid_y), abs(V) <= clip) %>%
    mutate(bin = floor(V / bin_width) * bin_width + bin_width / 2) %>%
    group_by(bin) %>%
    summarise(
      mean_y  = mean(resid_y),
      se_y    = sd(resid_y) / sqrt(n()),
      n       = n(),
      .groups = "drop"
    ) %>%
    filter(n >= min_n) %>%
    mutate(
      ci99_lo = mean_y - 2.576 * se_y,
      ci99_hi = mean_y + 2.576 * se_y
    )
}

plot_fig2 <- function(bin_df, ylab, title, y_clip = NULL) {
  left  <- bin_df %>% filter(bin <  0)
  right <- bin_df %>% filter(bin >= 0)
  
  p <- ggplot(bin_df, aes(x = bin, y = mean_y)) +
    # Smooth CI bands from the linear fit — one per side
    geom_smooth(data = left,  aes(x = bin, y = mean_y),
                method = "lm", formula = y ~ x,
                color = "black", fill = "grey70",
                se = TRUE, level = 0.99, linewidth = 1) +
    geom_smooth(data = right, aes(x = bin, y = mean_y),
                method = "lm", formula = y ~ x,
                color = "black", fill = "grey70",
                se = TRUE, level = 0.99, linewidth = 1) +
    geom_point(size = 1.5, color = "grey20") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_x_continuous(breaks = seq(-50, 50, by = 10)) +
    labs(x = "V: Profits Minus Losses", y = ylab, title = title) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor   = element_blank(),
          panel.background   = element_rect(fill = "grey93"),
          plot.background    = element_rect(fill = "white"))
  
  # Clip y-axis if requested, without dropping points from regression
  if (!is.null(y_clip))
    p <- p + coord_cartesian(ylim = c(-y_clip, y_clip))
  
  p
}

# Re-run everything
resid_refund_2002 <- residualize("TaxRefund",  df_2002, year_fe = TRUE)
resid_invest_2002 <- residualize("investment", df_2002, year_fe = TRUE)

bins_refund <- make_fig2_bins(resid_refund_2002)
bins_invest <- make_fig2_bins(resid_invest_2002)

# Re-draw with y-axis clipped to match paper scale
fig2_left  <- plot_fig2(bins_refund,
                        ylab   = "Tax Refund (Residual)",
                        title  = "Tax Refund (2002 Period)",
                        y_clip = 10)

fig2_right <- plot_fig2(bins_invest,
                        ylab   = "Investment (Residual)",
                        title  = "Investment (2002 Period)",
                        y_clip = 50)

# ggsave("figure2.png",
#        gridExtra::grid.arrange(fig2_left, fig2_right, ncol = 2,
#                                top = "Figure 2: Examples of Kinks in Variables"),
#        width = 12, height = 5, dpi = 150)
# print(fig2_left)
# print(fig2_right)

# fig2_left  <- plot_fig2(bins_refund,
#                         ylab  = "Tax Refund (Residual)",
#                         title = "Tax Refund (2002 Period)")
# fig2_right <- plot_fig2(bins_invest,
#                         ylab  = "Investment (Residual)",
#                         title = "Investment (2002 Period)")

# ggsave("figure2.png",
#        gridExtra::grid.arrange(fig2_left, fig2_right, ncol = 2,
#                                top = "Figure 2: Examples of Kinks in Variables"),
#        width = 12, height = 5, dpi = 150)
# print(fig2_left)
# print(fig2_right)

fig2 <- grid.arrange(fig2_left, fig2_right, ncol = 2,
                     top = "Figure 2: Examples of Kinks in Variables")

ggsave("figure2.png",
        grid.arrange(fig2_left, fig2_right, ncol = 2,
                     top = "Figure 2: Examples of Kinks in Variables"),
        width = 12, height = 5, dpi = 300)

#ggsave("figure2.png", fig2, width = 12, height = 5, dpi = 150)
print(fig2)

# # Combine side by side
# library(patchwork)
# fig2 <- fig2_left + fig2_right +
#   plot_annotation(title = "Figure 2: Examples of Kinks in Variables")
# 
# ggsave("figure2.png", fig2, width = 12, height = 5, dpi = 150)
# cat("Figure 2 saved.\n")

# ---- FIGURE 3: Histogram, $0.25M bins, clipped to [-15, 15] ----
make_fig3_hist <- function(df, title) {
  df %>%
    filter(!is.na(V), abs(V) <= 15) %>%
    ggplot(aes(x = V)) +
    geom_histogram(
      binwidth = 0.25,
      boundary = 0,
      fill     = "grey20",
      color    = "white",
      linewidth = 0.1
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    scale_x_continuous(breaks = seq(-15, 15, by = 5)) +
    labs(x = "V: Profits Minus Losses",
         y = "Frequency",
         title = title) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor  = element_blank(),
          panel.background  = element_rect(fill = "grey93"),
          plot.background   = element_rect(fill = "white"))
}

fig3_left  <- make_fig3_hist(df_2002, "2002 Policy")
#fig3_right <- make_fig3_hist(df_2009, "2009 Policy")

# fig3 <- grid.arrange(fig3_left, fig3_right, ncol = 2,
#                      top = "Figure 3: Distribution of Sample Firms around the Kink Point (V=0)")

ggsave("figure3.png",
       grid.arrange(fig3_left,
                    top = "Figure 3: Distribution of Sample Firms around the Kink Point (V=0)"),
       width = 7, height = 5, dpi = 300)

# ggsave("figure3.png",fig3_left,
#       width = 7, height = 5, dpi = 150)

print(fig3_left)

cat("Figures saved.\n")

# ---- 15. LATEX TABLES — write to .tex files --------------------

write_latex_table <- function(content, filename) {
  writeLines(content, filename)
  cat(sprintf("Written: %s\n", filename))
}

# ---- TABLE 4 ---------------------------------------------------
t4_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{First Stage: Average Firm Tax Rate}",
  "\\label{tab:table4}",
  "\\begin{tabular}{lcccc}",
  "\\hline\\hline",
  "Period & Coeff & SE & N & Paper \\\\",
  "\\hline"
)

t4_data <- list(
  list(period="2002 Policy", m=fs_2002, paper_coef="0.337", paper_sig="***"),
  list(period="2009 Policy", m=fs_2009, paper_coef="0.309", paper_sig="***")
)
for (row in t4_data) {
  p  <- tryCatch(pvalue(row$m)[["D_V"]], error=function(e) NA_real_)
  sg <- case_when(is.na(p)~"", p<.01~"***", p<.05~"**", p<.10~"*", TRUE~"")
  t4_lines <- c(t4_lines,
                sprintf("%s & %.3f%s & [%.3f] & %d & %s%s \\\\",
                        row$period,
                        coef(row$m)["D_V"], sg,
                        se(row$m)["D_V"],
                        row$m$nobs,
                        row$paper_coef, row$paper_sig))
}

t4_lines <- c(t4_lines,
              "\\hline\\hline",
              "\\multicolumn{5}{l}{\\footnotesize SEs in brackets, clustered by FF48 industry.} \\\\",
              "\\multicolumn{5}{l}{\\footnotesize *** $p<0.01$, ** $p<0.05$, * $p<0.10$} \\\\",
              "\\end{tabular}",
              "\\end{table}"
)
write_latex_table(t4_lines, "table4.tex")

# ---- TABLE 5 ---------------------------------------------------
# Exact values from paper Table 5
t5_paper_vals <- list(
  "Investment"   = list(coef= 0.403,   se=0.155,  sig="***"),
  "Δ Cash"       = list(coef=-0.189,   se=0.252,  sig=""),
  "Δ Total Debt" = list(coef=-0.606,   se=0.871,  sig=""),
  "Payout"       = list(coef=-0.0855,  se=0.130,  sig=""),
  "Other Uses"   = list(coef= 0.358,   se=0.579,  sig=""),
  "Δ Employment" = list(coef= 0.00271, se=0.00711,sig="")
)

t5_tex_labels <- c(
  "Investment"   = "Investment",
  "Δ Cash"       = "$\\Delta$ Cash",
  "Δ Total Debt" = "$\\Delta$ Total Debt",
  "Payout"       = "Payout",
  "Other Uses"   = "Other Uses",
  "Δ Employment" = "$\\Delta$ Employment"
)

t5_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Tax Refund Allocation: 2002 Policy Period (Refund Years 2002 and 2003)}",
  "\\label{tab:table5}",
  "\\begin{tabular}{lccclcc}",
  "\\hline\\hline",
  " & \\multicolumn{3}{c}{Replication} & & \\multicolumn{2}{c}{Paper} \\\\",
  "\\cline{2-4} \\cline{6-7}",
  "Outcome & Coeff & SE & N & & Coeff & SE \\\\",
  "\\hline"
)

for (i in seq_len(nrow(t5_rows))) {
  r       <- t5_rows[i, ]
  tex_lab <- t5_tex_labels[r$label]
  pv      <- t5_paper_vals[[r$label]]
  if (is.na(r$coef)) {
    t5_lines <- c(t5_lines,
                  sprintf("%s & \\multicolumn{3}{c}{[failed]} & & %.4f%s & [%.4f] \\\\",
                          tex_lab, pv$coef, pv$sig, pv$se))
  } else {
    t5_lines <- c(t5_lines,
                  sprintf("%s & %.3f%s & [%.3f] & %d & & %.4f%s & [%.4f] \\\\",
                          tex_lab, r$coef, r$sig, r$se, r$n,
                          pv$coef, pv$sig, pv$se))
  }
}

t5_lines <- c(t5_lines,
              "\\hline\\hline",
              "\\multicolumn{7}{l}{\\footnotesize SEs in brackets, clustered by FF48 industry.} \\\\",
              "\\multicolumn{7}{l}{\\footnotesize *** $p<0.01$, ** $p<0.05$, * $p<0.10$} \\\\",
              "\\end{tabular}",
              "\\end{table}"
)
write_latex_table(t5_lines, "table5.tex")

# ---- TABLE 6 ---------------------------------------------------
# Exact values from paper Table 6
t6_paper_vals <- list(
  # Panel A
  "KZ: constrained (top Q)"    = list(coef=1.091,  se=0.659, sig="*"),
  "KZ: unconstrained (bot Q)"  = list(coef=0.0534, se=0.228, sig=""),
  "Payout: zero (constrained)" = list(coef=1.004,  se=0.535, sig="*"),
  "Payout: >0 (unconstrained)" = list(coef=0.418,  se=0.307, sig=""),
  "DFF: ≤0 (constrained)"      = list(coef=1.002,  se=0.509, sig="**"),
  "DFF: >0 (unconstrained)"    = list(coef=0.501,  se=0.389, sig=""),
  # Panel B
  "High Tobin's Q (top Q)"     = list(coef= 1.017,  se=0.593, sig="*"),
  "Low Tobin's Q (bot Q)"      = list(coef=-1.454,  se=1.481, sig=""),
  "Refund year: 2002"          = list(coef= 0.277,  se=0.384, sig=""),
  "Refund year: 2003"          = list(coef= 0.753,  se=0.243, sig="***"),
  "Multinational"              = list(coef= 0.676,  se=0.372, sig="*"),
  "Domestic"                   = list(coef= 0.171,  se=0.177, sig="")
)

t6_tex_labels <- c(
  "KZ: constrained (top Q)"    = "KZ: constrained (top quartile)",
  "KZ: unconstrained (bot Q)"  = "KZ: unconstrained (bottom quartile)",
  "Payout: zero (constrained)" = "Payout: zero (constrained)",
  "Payout: >0 (unconstrained)" = "Payout: $>$0 (unconstrained)",
  "DFF: ≤0 (constrained)"      = "DFF: $\\leq$0 (constrained)",
  "DFF: >0 (unconstrained)"    = "DFF: $>$0 (unconstrained)",
  "High Tobin's Q (top Q)"     = "High Tobin's Q (top quartile)",
  "Low Tobin's Q (bot Q)"      = "Low Tobin's Q (bottom quartile)",
  "Refund year: 2002"          = "Refund year: 2002",
  "Refund year: 2003"          = "Refund year: 2003",
  "Multinational"              = "Multinational",
  "Domestic"                   = "Domestic"
)

build_t6_panel <- function(rows, panel_title) {
  out <- c(
    sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", panel_title),
    "\\hline"
  )
  for (i in seq_len(nrow(rows))) {
    r       <- rows[i, ]
    tex_lab <- t6_tex_labels[r$label]
    pv      <- t6_paper_vals[[r$label]]
    if (is.na(r$coef)) {
      out <- c(out,
               sprintf("%s & \\multicolumn{3}{c}{[failed]} & & %.3f%s & [%.3f] \\\\",
                       tex_lab, pv$coef, pv$sig, pv$se))
    } else {
      out <- c(out,
               sprintf("%s & %.3f%s & [%.3f] & %d & & %.3f%s & [%.3f] \\\\",
                       tex_lab, r$coef, r$sig, r$se, r$n,
                       pv$coef, pv$sig, pv$se))
    }
  }
  out
}

t6_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Investment by Subgroup: 2002 Policy Period}",
  "\\label{tab:table6}",
  "\\begin{tabular}{lccclcc}",
  "\\hline\\hline",
  " & \\multicolumn{3}{c}{Replication} & & \\multicolumn{2}{c}{Paper} \\\\",
  "\\cline{2-4} \\cline{6-7}",
  "Subgroup & Coeff & SE & N & & Coeff & SE \\\\",
  "\\hline",
  build_t6_panel(rows_a, "Panel A: Financial Constraints"),
  "\\hline",
  build_t6_panel(rows_b, "Panel B: Investment Opportunities"),
  "\\hline\\hline",
  "\\multicolumn{7}{l}{\\footnotesize SEs in brackets, clustered by FF48 industry.} \\\\",
  "\\multicolumn{7}{l}{\\footnotesize *** $p<0.01$, ** $p<0.05$, * $p<0.10$} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)
write_latex_table(t6_lines, "table6.tex")
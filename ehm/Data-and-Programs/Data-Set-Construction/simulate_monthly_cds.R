##############################################################################
#### Simulate monthly CDS panel when original "monthly cds.dta" is missing ####
##############################################################################

suppressPackageStartupMessages({
  library(haven)
  library(readxl)
})

set.seed(12345)

# Targets from user
TARGET_MEAN <- 0.0083  # 0.83%
TARGET_SD   <- 0.0088  # 0.88%
TARGET_MIN  <- 0.0005  # 0.05%
TARGET_MAX  <- 0.0547  # 5.47%

start_date <- as.Date("2002-01-01")
end_date   <- as.Date("2013-12-01")
months_seq <- seq.Date(start_date, end_date, by = "month")
n_t <- length(months_seq)

# Use shortnames that match create_deposit_data_final.do and demand merges.
banks <- data.frame(
  shortname = c(
    "BofA", "JPM", "Wells", "US", "Citi", "PNC", "Wachovia", "RBS",
    "Suntrust", "TD", "Fifth", "Regions", "HSBC", "BB&T", "M&T",
    "NatCity", "Key", "Santander"
  ),
  certno = c(
    3510, 628, 3511, 6548, 7213, 6384, 33869, 57957,
    867, 18409, 6672, 12368, 57890, 9846, 588, 6557, 17534, 29950
  ),
  stringsAsFactors = FALSE
)
n_b <- nrow(banks)

to_stata_month <- function(d) {
  lt <- as.POSIXlt(d)
  (lt$year + 1900 - 1960) * 12 + lt$mon
}

parse_reference <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- tryCatch(read_excel(path), error = function(e) NULL)
  if (is.null(x) || ncol(x) < 2) return(NULL)

  # Guess date column
  date_col <- NULL
  for (j in seq_len(ncol(x))) {
    v <- x[[j]]
    if (inherits(v, "Date") || inherits(v, "POSIXct")) {
      date_col <- j
      break
    }
    parsed <- suppressWarnings(as.Date(v))
    if (mean(!is.na(parsed)) > 0.7) {
      date_col <- j
      break
    }
  }
  if (is.null(date_col)) date_col <- 1

  # Guess CDS column (numeric, most variation)
  num_idx <- which(vapply(x, is.numeric, logical(1)))
  num_idx <- setdiff(num_idx, date_col)
  if (length(num_idx) == 0) return(NULL)
  sds <- vapply(num_idx, function(j) sd(x[[j]], na.rm = TRUE), numeric(1))
  cds_col <- num_idx[which.max(sds)]

  d <- x[[date_col]]
  if (!inherits(d, "Date")) d <- suppressWarnings(as.Date(d))
  v <- as.numeric(x[[cds_col]])
  keep <- !is.na(d) & !is.na(v)
  d <- d[keep]
  v <- v[keep]
  if (length(v) < 12) return(NULL)

  # Convert scale to decimal rate if needed
  q99 <- as.numeric(quantile(v, 0.99, na.rm = TRUE))
  if (q99 > 20) {
    v <- v / 10000  # bps -> decimal
  } else if (q99 > 1.5) {
    v <- v / 100    # percent points -> decimal
  }

  data.frame(date = as.Date(format(d, "%Y-%m-01")), spread = v)
}

calc_amp <- function(ref_df, fallback_core = 0.012, fallback_dist = 0.028) {
  if (is.null(ref_df) || nrow(ref_df) == 0) return(c(core = fallback_core, dist = fallback_dist))
  pre <- ref_df$spread[ref_df$date < as.Date("2008-01-01")]
  crisis <- ref_df$spread[ref_df$date >= as.Date("2008-01-01") & ref_df$date <= as.Date("2010-12-01")]
  if (length(pre) < 3 || length(crisis) < 3) return(c(core = fallback_core, dist = fallback_dist))
  c(core = max(0.003, mean(crisis, na.rm = TRUE) - mean(pre, na.rm = TRUE)),
    dist = max(0.010, max(crisis, na.rm = TRUE) - mean(pre, na.rm = TRUE)))
}

# Reference files (used for crisis intensity calibration if readable)
jpm_ref  <- parse_reference("data/jpm08.xlsx")
bofa_ref <- parse_reference("data/bofa08.xlsx")
citi_ref <- parse_reference("data/citi08.xlsx")

core_amp <- mean(c(calc_amp(jpm_ref)["core"], calc_amp(bofa_ref)["core"]))
dist_amp <- calc_amp(citi_ref)["dist"]

# Common time factor with persistence and crisis pulses
eps_t <- rnorm(n_t, mean = 0, sd = 0.0008)
common <- numeric(n_t)
for (t in 2:n_t) common[t] <- 0.92 * common[t - 1] + eps_t[t]

is_pre2008 <- months_seq < as.Date("2008-01-01")
pre_shift <- ifelse(is_pre2008, -0.0016, 0.0000)

g1 <- exp(-0.5 * ((as.numeric(months_seq - as.Date("2008-10-01")) / 30.5) / 4.5)^2)
g2 <- exp(-0.5 * ((as.numeric(months_seq - as.Date("2009-03-01")) / 30.5) / 6.5)^2)
crisis_shape <- 0.65 * g1 + 0.35 * g2
crisis_shape <- crisis_shape / max(crisis_shape)

base_level <- c(
  BofA = 0.0048, JPM = 0.0041, Wells = 0.0037, US = 0.0032, Citi = 0.0056,
  PNC = 0.0038, Wachovia = 0.0048, RBS = 0.0050, Suntrust = 0.0044,
  TD = 0.0032, Fifth = 0.0048, Regions = 0.0047, HSBC = 0.0041,
  `BB&T` = 0.0037, `M&T` = 0.0034, NatCity = 0.0049, Key = 0.0042, Santander = 0.0045
)

# Distress multipliers: stronger spikes for known stressed banks.
distress_mult <- c(
  BofA = 0.55, JPM = 0.40, Wells = 0.48, US = 0.35, Citi = 1.00,
  PNC = 0.45, Wachovia = 0.95, RBS = 0.85, Suntrust = 0.62,
  TD = 0.35, Fifth = 0.68, Regions = 0.70, HSBC = 0.60,
  `BB&T` = 0.50, `M&T` = 0.40, NatCity = 0.88, Key = 0.62, Santander = 0.52
)

idiosyncratic <- matrix(0, nrow = n_t, ncol = n_b)
for (b in seq_len(n_b)) {
  sig_b <- 0.0007 + 0.0006 * distress_mult[banks$shortname[b]]
  e_b <- rnorm(n_t, mean = 0, sd = sig_b)
  for (t in 2:n_t) idiosyncratic[t, b] <- 0.86 * idiosyncratic[t - 1, b] + e_b[t]
}

sim_panel <- vector("list", n_b)
for (b in seq_len(n_b)) {
  s <- banks$shortname[b]
  level <- base_level[s]
  crisis <- (core_amp + (dist_amp - core_amp) * distress_mult[s]) * crisis_shape
  raw <- level + pre_shift + common + idiosyncratic[, b] + crisis
  sim_panel[[b]] <- data.frame(
    certno = banks$certno[b],
    shortname = s,
    month = to_stata_month(months_seq),
    spread5y_raw = raw,
    stringsAsFactors = FALSE
  )
}

cds <- do.call(rbind, sim_panel)

# Match global moments with a linear transform + clipping.
match_moments <- function(x, m_t, s_t, lo, hi, iter = 5) {
  y <- x
  for (k in seq_len(iter)) {
    m <- mean(y, na.rm = TRUE)
    s <- sd(y, na.rm = TRUE)
    if (is.na(s) || s == 0) break
    b <- s_t / s
    a <- m_t - b * m
    y <- a + b * y
    y <- pmin(hi, pmax(lo, y))
  }
  y
}

cds$spread5y <- match_moments(cds$spread5y_raw, TARGET_MEAN, TARGET_SD, TARGET_MIN, TARGET_MAX)

# Ensure at least one observation hits the requested max and min.
peak_idx <- which.max(cds$spread5y_raw + 10 * (cds$shortname == "Citi"))
trough_idx <- which.min(cds$spread5y_raw)
cds$spread5y[peak_idx] <- TARGET_MAX
cds$spread5y[trough_idx] <- TARGET_MIN

# Create 1Y spread as a noisy affine transform of 5Y.
noise_1y <- rnorm(nrow(cds), 0, 0.00045)
cds$spread1y <- pmin(TARGET_MAX, pmax(TARGET_MIN, 0.87 * cds$spread5y + 0.00025 + noise_1y))

cds <- cds[order(cds$shortname, cds$month), c("certno", "shortname", "month", "spread1y", "spread5y")]
row.names(cds) <- NULL

out_dta <- "Data-and-Programs/Data-Sets/monthly cds.dta"
out_csv <- "Data-and-Programs/Data-Sets/monthly_cds_simulated.csv"

write_dta(cds, out_dta)
write.csv(cds, out_csv, row.names = FALSE)

cat("\nSaved simulated monthly CDS:\n")
cat("  ", out_dta, "\n")
cat("  ", out_csv, "\n")
cat("\nSummary (spread5y):\n")
print(summary(cds$spread5y))
cat("\nMean:", mean(cds$spread5y), " SD:", sd(cds$spread5y),
    " Min:", min(cds$spread5y), " Max:", max(cds$spread5y), "\n")

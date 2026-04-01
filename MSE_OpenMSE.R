# =============================================================
# MSE Antarctic Krill — openMSE / MSEtool
# HCR 1: Constant Catch 145,000 t (CCAMLR current level)
# HCR 2: Proportional Catch to observed biomass (5% B_obs)
# OM Grid: sigmaR × R0_shift (climate scenarios)
# Version: 7.0 — maxage=8 biological, CC145, Find+HCR bugs fixed
# =============================================================
#
# ═══════════════════════════════════════════════════════════
# CORRECTIONS APPLIED IN v7.0
# ═══════════════════════════════════════════════════════════
#
# FIX 1 — maxage=8 biological (real krill, no phantom ages)
#   Age arrays: [nsim, 9, total_yr]  (ages 0:8 → maxage+1 cols)
#   Perr_y:     [nsim, 66]            (8 + 28 + 30)
#
# FIX 2 — Find: pull(`F:_Bio_smry`) does not exist in SS3
#   Real columns: F:_1 ... F:_5  → grep + rowSums
#   Fallback: catch/Bio_smry if F:_ columns also missing
#
# FIX 3 — HCR_proportional: Data@Ind is a relative index (~0-2)
#   B_abs = B_rel * B0  → TAC = 0.05 * B_abs (correct scale)
#
# ═══════════════════════════════════════════════════════════

library(MSEtool)
library(r4ss)
library(tidyverse)
library(patchwork)

# -------------------------------------------------------------
# GLOBAL PARAMETERS
# -------------------------------------------------------------

ss3_dir     <- "s1.4"
yr_start    <- 1991
yr_final    <- 2018
n_sim       <- 48
n_proyear   <- 30
sigmaR_base <- 1.2
AC_base     <- 0.5

# Krill biological parameters
maxage  <- 8        # true biological maximum age of krill
M_base  <- 0.27     # natural mortality adults (yr-1)
h       <- 0.85     # Beverton-Holt steepness
R0      <- 368857000
B0      <- 2028880
D_2018  <- 0.8085

# Von Bertalanffy
Linf    <- 60.0
K       <- 0.48
t0      <- -0.1

# Length-weight  W = a*L^b  (W in tonnes, L in mm)
a_LW    <- 2.36e-12
b_LW    <- 3.19

# Maturity and selectivity
L50_mat <- 35.0; L95_mat <- 45.0
L50_sel <- 30.0; L95_sel <- 45.0

nyears   <- yr_final - yr_start + 1   # 28
total_yr <- nyears + n_proyear         # 58

# ── Array dimensions for MSEtool ────────────────────────────
# MSEtool indexes ages 0:maxage → needs maxage+1 columns
# Perr_y needs [nsim, maxage + nyears + proyears]

n_ages_mse <- maxage + 1                          # = 9 (ages 0:8)
perr_ncol  <- maxage + nyears + n_proyear          # = 66

cat("=== Array configuration ===\n")
cat("Biological maxage:         ", maxage, "\n")
cat("n_ages_mse (age cols):     ", n_ages_mse, "(ages 0 to", maxage, ")\n")
cat("Perr_y ncols:              ", perr_ncol,
    "(=", maxage, "+", nyears, "+", n_proyear, ")\n")
cat("Age arrays:                [", n_sim, "×", n_ages_mse, "×", total_yr, "]\n\n")

# -------------------------------------------------------------
# STEP 1: READ SS3
# -------------------------------------------------------------

rep <- SS_output(dir = ss3_dir, verbose = FALSE, printstats = FALSE)

# Historical SSB/B0 time series from SS3 (reference line in wormplot)
hist_ssb <- rep$timeseries |>
  filter(Yr >= yr_start, Yr <= yr_final) |>
  arrange(Yr) |>
  transmute(year = Yr, ssb_rel = SpawnBio / B0)

rec_devs <- rep$recruit |>
  filter(Yr >= yr_start, Yr <= yr_final) |>
  arrange(Yr) |>
  pull(dev)

rec_devs[is.na(rec_devs)] <- 0

cat("Recruitment deviations:", length(rec_devs), "years\n")
cat("Dev range:", round(range(rec_devs), 3), "\n")

# Historical Perr: multipliers (linear scale, not log)
Perr_hist_vec <- exp(rec_devs)

# Historical recruitment from SS3 — use exp_recr (expected recruits in numbers)
rec_numbers_ss3 <- rep$recruit |>
  filter(Yr >= yr_start, Yr <= yr_final) |>
  arrange(Yr) |>
  pull(exp_recr)  # Expected recruits (individuals) from SS3

# Convert to biomass of age-0: numbers × weight-at-age-0
rec_biomass_ss3 <- rec_numbers_ss3 * W_all[1]  # tonnes

# Data frame for historical recruitment reference line
hist_rec <- data.frame(
  year = seq(yr_start, yr_final),
  rec_tonnes = rec_biomass_ss3 / 1e3  # convert to thousands of tonnes
)

cat("Historical recruitment (age-0 biomass, kt):\n")
cat("  Mean:", round(mean(hist_rec$rec_tonnes), 0), "kt\n")
cat("  Range:", round(min(hist_rec$rec_tonnes), 0), "–", round(max(hist_rec$rec_tonnes), 0), "kt\n")

# Historical matrix: nsim × nyears — small inter-sim variability
Perr_hist <- matrix(rep(Perr_hist_vec, n_sim),
                    nrow = n_sim, ncol = nyears, byrow = TRUE)
set.seed(123)
noise <- matrix(rnorm(n_sim * nyears, 0, 0.05), nrow = n_sim)
Perr_hist <- Perr_hist * exp(noise)
cat("Perr_hist dim:", paste(dim(Perr_hist), collapse = "×"), "\n")

# -------------------------------------------------------------
# STEP 2: BIOLOGICAL ARRAYS (ages 0:maxage = 0:8)
# -------------------------------------------------------------
# FIX 1: ages 0 to maxage → n_ages_mse = 9 columns (maxage+1)

ages_all <- 0:maxage   # 0, 1, 2, ..., 8

# Von Bertalanffy for all ages
L_all <- Linf * (1 - exp(-K * (ages_all - t0)))

# Weight (tonnes)
W_all <- a_LW * L_all^b_LW

# Maturity
Mat_all <- 1 / (1 + exp(-log(19) * (L_all - L50_mat) / (L95_mat - L50_mat)))

# Selectivity
Sel_all <- 1 / (1 + exp(-log(19) * (L_all - L50_sel) / (L95_sel - L50_sel)))

# Age-specific natural mortality (9 values for ages 0:8)
# Standard krill biological pattern: high post-larval mortality, decreasing with age
M_all <- c(
  M_base * 3.0,   # age 0 — post-larva: very high M
  M_base * 2.5,   # age 1 — early juvenile
  M_base * 2.0,   # age 2 — late juvenile
  M_base * 1.5,   # age 3 — sub-adult
  M_base * 1.2,   # age 4 — pre-adult
  M_base * 1.0,   # age 5 — adult
  M_base * 1.0,   # age 6 — adult
  M_base * 1.0,   # age 7 — adult
  M_base * 1.0    # age 8 — adult (maxage)
)
stopifnot(length(M_all) == n_ages_mse)   # must be 9

cat("\nBiology by age (ages 0:", maxage, "):\n")
df_bio <- data.frame(
  age  = ages_all,
  L_mm = round(L_all, 1),
  W_g  = round(W_all * 1e9, 4),   # in micrograms (approx)
  Mat  = round(Mat_all, 3),
  Sel  = round(Sel_all, 3),
  M    = round(M_all, 3)
)
print(df_bio)

# Function to build 3D array [nsim, n_ages_mse, total_yr]
make_age_array <- function(vec, n_sim, total_yr) {
  arr <- array(NA, dim = c(n_sim, length(vec), total_yr))
  for (s in 1:n_sim)
    for (t in 1:total_yr)
      arr[s, , t] <- vec
  arr
}

Wt_age_arr  <- make_age_array(W_all,   n_sim, total_yr)
Mat_age_arr <- make_age_array(Mat_all, n_sim, total_yr)
M_age_arr   <- make_age_array(M_all,   n_sim, total_yr)
V_arr       <- make_age_array(Sel_all, n_sim, total_yr)

cat("\nAge arrays — dim:", paste(dim(Wt_age_arr), collapse = "×"),
    "= [nsim × (maxage+1) × total_yr] ✓\n")

# ── R0 adjusted for consistency with SS3 B0 ──────────────────────────────
# PROBLEM: MSEtool computes B0_internal = R0 × virgin_SPR
#   Using R0_ss3 directly → B0_internal ≈ 79 t (not 2M t) → immediate collapse
# SOLUTION: R0_mse = B0_ss3 / SPR_0  →  MSEtool sees the correct B0
#
# SPR with plus group at maxage (MSEtool internal convention)
l_spr <- c(1, cumprod(exp(-M_all[-length(M_all)])))
l_spr[length(l_spr)] <- l_spr[length(l_spr)] / (1 - exp(-M_all[length(M_all)]))
SPR_0  <- sum(l_spr * W_all * Mat_all)   # tonnes per virgin recruit
R0_mse <- round(B0 / SPR_0)              # R0 consistent with B0_ss3

cat(sprintf("\nVirgin SPR (plus group): %.4e t/recruit\n", SPR_0))
cat(sprintf("R0 SS3:      %.4e ind  → B0_internal MSEtool = %.1f t  [WRONG]\n",
            R0, R0 * SPR_0))
cat(sprintf("Adjusted R0: %.4e ind  → B0_internal MSEtool = %.0f t  [= B0_ss3 ✓]\n",
            R0_mse, R0_mse * SPR_0))

# -------------------------------------------------------------
# STEP 3: sim_perr_proj (AR1 log-normal)
# -------------------------------------------------------------

sim_perr_proj <- function(n_sim, n_proy, sigmaR, AC, seed = 42) {
  set.seed(seed)
  mat <- matrix(NA, nrow = n_sim, ncol = n_proy)
  for (s in seq_len(n_sim)) {
    eps      <- numeric(n_proy)
    eps_prev <- 0
    for (t in seq_len(n_proy)) {
      eps[t]   <- AC * eps_prev + sqrt(1 - AC^2) * rnorm(1, 0, sigmaR)
      eps_prev <- eps[t]
    }
    mat[s, ] <- exp(eps)
  }
  mat
}

# -------------------------------------------------------------
# STEP 4: Perr_y with correct dimensions [nsim, 66]
# -------------------------------------------------------------
# FIX 1: Perr_y = [nsim, maxage + nyears + proyears]
#         = [48, 8 + 28 + 30] = [48, 66]
#
# Block 1: pre-history (cols 1:maxage)          → cohorts pre-1991 at start
# Block 2: historical  (cols maxage+1 to +nyears) → SS3 deviations
# Block 3: projection  (remaining cols)          → AR1 log-normal

build_perr_y <- function(sigmaR_i, AC_i, seed_i = 42) {
  # Pre-history: maxage columns (~1 ± small noise)
  set.seed(seed_i + 1000)
  Perr_pre  <- matrix(
    exp(rnorm(n_sim * maxage, mean = 0, sd = 0.1)),
    nrow = n_sim, ncol = maxage    # [48 × 8]
  )
  # Projection: AR1 with sigmaR and AC from scenario
  Perr_proj <- sim_perr_proj(n_sim, n_proyear, sigmaR_i, AC_i, seed = seed_i)
  # Combine: [48, 8] | [48, 28] | [48, 30] = [48, 66]
  cbind(Perr_pre, Perr_hist, Perr_proj)
}

Perr_y_base <- build_perr_y(sigmaR_base, AC_base)
cat("\nPerr_y_base dim:", paste(dim(Perr_y_base), collapse = "×"),
    "(expected: 48×", perr_ncol, ")\n")
stopifnot(all(dim(Perr_y_base) == c(n_sim, perr_ncol)))

# -------------------------------------------------------------
# STEP 5: EXTRACT HISTORICAL Find (FIX 2)
# -------------------------------------------------------------
# Real columns in Report.sso: F:_1, F:_2, ..., F:_5
# pull(`F:_Bio_smry`) DOES NOT EXIST → dplyr error

ts_filt <- rep$timeseries |>
  filter(Yr >= yr_start, Yr <= yr_final) |>
  arrange(Yr)

# Search for fleet-specific F columns: pattern "^F:_"
f_fleet_cols <- grep("^F:_", names(ts_filt), value = TRUE)

if (length(f_fleet_cols) > 0) {
  # Sum of F across all fleets = total F
  Find_ss3 <- rowSums(ts_filt[, f_fleet_cols, drop = FALSE], na.rm = TRUE)
  cat("\nFind extracted from", length(f_fleet_cols), "fleet columns:",
      paste(f_fleet_cols, collapse = ", "), "\n")
} else {
  # Fallback: F ≈ Catch / Biomass (if F:_ columns not found)
  cat("⚠ No F:_ columns found — using fallback catch/biomass\n")
  obs_cols    <- grep("^obs_cat:_", names(ts_filt), value = TRUE)
  total_catch <- if (length(obs_cols) > 0)
    rowSums(ts_filt[, obs_cols, drop = FALSE], na.rm = TRUE)
  else ts_filt$`obs_cat:_1`
  Find_ss3 <- pmin(total_catch / pmax(ts_filt$Bio_smry, 1), 2.0)
}

Find_ss3[is.na(Find_ss3)] <- 0
Find_mat <- matrix(rep(Find_ss3, n_sim), nrow = n_sim, byrow = TRUE)

cat("Historical Find — range:", round(range(Find_mat), 4), "\n")
cat("Historical Find — mean:", round(mean(Find_mat), 4), "\n")

# -------------------------------------------------------------
# STEP 6: BUILD BASE OM
# -------------------------------------------------------------

# Use testOM as template — all required slots are already valid.
# We only overwrite the krill-specific parameters.
data("testOM")
om <- testOM

# ── Basic configuration ─────────────────────────────────────
om@nsim      <- n_sim
om@nyears    <- nyears
om@proyears  <- n_proyear
om@maxage    <- maxage
om@CurrentYr <- yr_final
om@interval  <- 1L
om@SRrel     <- 1L   # Beverton-Holt

# ── Biology (overwrites testOM values) ─────────────────────
om@R0    <- R0_mse   # adjusted: B0_ss3 / SPR_0 (not R0_ss3 directly)
om@M     <- c(M_base, M_base)
om@Msd   <- c(0, 0)
om@h     <- c(h, h)
om@D     <- c(D_2018, D_2018)
om@AC    <- c(AC_base, AC_base)
om@Perr  <- c(sigmaR_base, sigmaR_base)

# Von Bertalanffy
om@Linf   <- c(Linf, Linf)
om@Linfsd <- c(0, 0)
om@K      <- c(K, K)
om@Ksd    <- c(0, 0)
om@t0     <- c(t0, t0)
om@LenCV  <- c(0.1, 0.1)

# Length-weight — single scalar (MSEtool requires length 1)
om@a <- a_LW
om@b <- b_LW

# Maturity
om@L50    <- c(L50_mat, L50_mat)
om@L50_95 <- c(L95_mat - L50_mat, L95_mat - L50_mat)

# Historical effort adjusted to model period
om@EffYears <- c(yr_start, yr_final)
om@EffLower <- c(1, 1)
om@EffUpper <- c(1, 1)

# Observation
om@Cobs      <- c(0.1, 0.1)
om@Iobs      <- c(0.2, 0.2)

# Implementation without error
om@TACFrac <- c(1, 1)
om@TACSD   <- c(0, 0)

# ── cpars: controls all dynamics via SS3 ───────────────────
om@cpars <- list(
  Wt_age     = Wt_age_arr,    # [48, 9, 58] — ages 0:8
  Mat_age    = Mat_age_arr,   # [48, 9, 58]
  M_ageArray = M_age_arr,     # [48, 9, 58]
  V          = V_arr,         # [48, 9, 58]
  Perr_y     = Perr_y_base,   # [48, 66]
  Find       = Find_mat       # [48, 28]
)

cat("\n=== Base OM constructed ===\n")
cat("om@maxage:", om@maxage, "\n")
cat("Wt_age dim:", paste(dim(om@cpars$Wt_age), collapse = "×"),
    "(expected:", paste(c(n_sim, n_ages_mse, total_yr), collapse = "×"), ")\n")
cat("Perr_y dim:", paste(dim(om@cpars$Perr_y), collapse = "×"),
    "(expected:", paste(c(n_sim, perr_ncol), collapse = "×"), ")\n")

# ── CheckOM ─────────────────────────────────────────────────
cat("\nRunning CheckOM...\n")
tryCatch({
  chk <- CheckOM(om)
  cat("CheckOM output:\n")
  print(chk)
}, error = function(e) {
  cat("CheckOM error:", conditionMessage(e), "\n")
})

# ── Quick OM base test ──────────────────────────────────────
cat("\nTesting base OM with runMSE('FMSYref')...\n")
test_om <- tryCatch(
  runMSE(OM = om, MPs = "FMSYref", parallel = FALSE),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(test_om)) {
  cat("✅ Base OM OK\n")
  cat("Median SSB year  1: ", round(median(test_om@SSB[, 1,  1])), "t\n")
  cat("Median SSB year 15:", round(median(test_om@SSB[, 1, 15])), "t\n")
  cat("Median SSB year 30:", round(median(test_om@SSB[, 1, 30])), "t\n")
  cat("Reference B0:       ", B0, "t\n")
  cat("Expected initial D: ", D_2018, " | Observed:",
      round(median(test_om@SSB[, 1, 1]) / B0, 3), "\n")
} else {
  cat("❌ Base OM failed. Review error messages above.\n")
  stop("Base OM not valid — stopping execution")
}

# -------------------------------------------------------------
# STEP 7: 9-OM GRID
# -------------------------------------------------------------

om_grid <- expand.grid(
  sigmaR   = c(0.8, 1.2, 1.5),
  R0_shift = c(0.0, -0.2, -0.4)
) |> mutate(om_id = paste0("sR", sigmaR, "_R0", R0_shift * 100))

cat("\n=== OM Grid ===\n")
print(om_grid)

om_list <- vector("list", nrow(om_grid))

for (i in seq_len(nrow(om_grid))) {
  om_i        <- om
  om_i@Perr   <- rep(om_grid$sigmaR[i], 2)
  om_i@R0     <- R0_mse * (1 + om_grid$R0_shift[i])  # scale adjusted R0
  om_i@cpars$Perr_y <- build_perr_y(om_grid$sigmaR[i], AC_base, seed_i = 42 + i)
  om_list[[i]] <- om_i
  cat("OM", i, "—", om_grid$om_id[i],
      "| sigmaR:", om_i@Perr[1],
      "| R0:", round(om_i@R0 / 1e6, 1), "M ind\n")
}

# -------------------------------------------------------------
# STEP 8: DEFINE MPs
# -------------------------------------------------------------

# ── HCR 1: Constant Catch 145,000 t ───────────────────────
CC145 <- function(x, Data, reps = 1, ...) {
  Rec     <- new("Rec")
  Rec@TAC <- rep(145000, reps)
  return(Rec)
}
class(CC145) <- "MP"

# ── HCR 2: Proportional to biomass (FIX 3) ────────────────
# Data@Ind is a relative index (~0-2), NOT biomass in tonnes
# B_abs = B_rel * B0  → TAC = 5% of B_abs
B0_MSE <- B0   # capture B0 before closure

HCR_proportional <- function(x, Data, reps = 1, ...) {
  ind   <- Data@Ind[x, ]
  B_rel <- ind[max(which(!is.na(ind)))]   # relative index (last non-NA value)
  B_abs <- B_rel * B0_MSE                  # scale to tonnes using B0
  TAC   <- max(B_abs * 0.05, 0)            # 5% estimated absolute biomass
  Rec   <- new("Rec")
  Rec@TAC <- rep(TAC, reps)
  return(Rec)
}
class(HCR_proportional) <- "MP"

MPs_names <- c("CC145", "HCR_proportional")
cat("\nMPs defined:", paste(MPs_names, collapse = ", "), "\n")
cat("CC145: Fixed TAC = 145,000 t\n")
cat("HCR_proportional: TAC = 5% × (B_rel × B0), B0 =", B0, "t\n")

# ══════════════════════════════════════════════════════════════
# PRE-MSE DIAGNOSTICS
# ══════════════════════════════════════════════════════════════
cat("\n====== PRE-MSE DIAGNOSTICS ======\n")

# ── 1. B0 consistency ───────────────────────────────────────
l_check <- c(1, cumprod(exp(-M_all[-length(M_all)])))
l_check[length(l_check)] <- l_check[length(l_check)] / (1 - exp(-M_all[length(M_all)]))
B0_check  <- R0_mse * sum(l_check * W_all * Mat_all)
SSB_init  <- D_2018 * B0_check

cat(sprintf("B0 SS3:            %10.0f t\n", B0))
cat(sprintf("B0 MSEtool:        %10.0f t  (should be ≈ B0 SS3)\n", B0_check))
cat(sprintf("Initial SSB:       %10.0f t  (D_2018 × B0)\n", SSB_init))
cat(sprintf("CC145/SSB_init:         %.4f  (should be << 1)\n", 145000/SSB_init))
cat(sprintf("Implicit F CC145:       %.4f  (adult M = %.2f)\n\n", 145000/SSB_init, M_base))

# ── 2. Virgin equilibrium table ─────────────────────────────
N_eq    <- R0_mse * l_check
SSB_age <- N_eq * W_all * Mat_all
df_eq   <- data.frame(
  age    = 0:maxage,
  l_a    = round(l_check, 4),
  N_bill = round(N_eq / 1e9, 1),
  W_g    = round(W_all * 1e6, 3),
  Mat    = round(Mat_all, 3),
  SSB_t  = round(SSB_age)
)
print(df_eq)
cat(sprintf("Total virgin SSB: %s t\n\n",
            format(round(sum(SSB_age)), big.mark = ",")))

# ── 3. Quick test: 3 sims, 5 years ─────────────────────────
n_proy_test  <- 5
n_tot_test   <- nyears + n_proy_test          # 33
perr_test    <- maxage + nyears + n_proy_test  # 41

om_test              <- om
om_test@nsim         <- 3
om_test@proyears     <- n_proy_test
om_test@cpars$Wt_age     <- om@cpars$Wt_age[1:3,     , 1:n_tot_test]
om_test@cpars$Mat_age    <- om@cpars$Mat_age[1:3,    , 1:n_tot_test]
om_test@cpars$M_ageArray <- om@cpars$M_ageArray[1:3, , 1:n_tot_test]
om_test@cpars$V          <- om@cpars$V[1:3,          , 1:n_tot_test]
om_test@cpars$Find       <- om@cpars$Find[1:3, ]
om_test@cpars$Perr_y     <- om@cpars$Perr_y[1:3, 1:perr_test]

res_test <- tryCatch(
  runMSE(OM = om_test, MPs = "FMSYref", parallel = FALSE),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(res_test)) {
  ssb_rel <- res_test@SSB[, 1, ] / B0_check
  cat(sprintf("SSB/B0 — year 1: %.3f  year 3: %.3f  year 5: %.3f\n",
              median(ssb_rel[, 1]), median(ssb_rel[, 3]), median(ssb_rel[, 5])))
  cat(sprintf("Expected t=1: ~%.3f (= D_2018)\n", D_2018))
  if (median(ssb_rel[, 1]) < 0.05)
    cat("❌ SSB/B0 ≈ 0 — scale problem, do not run full MSE\n") else
      cat("✅ OM validated — ready for MSE 9×2\n")
}
cat("══════════════════════════════════════════\n\n")

# -------------------------------------------------------------
# STEP 9: RUN MSE
# -------------------------------------------------------------

cat("\n=== Starting MSE — 9 OMs × 2 MPs ===\n")
mse_list <- vector("list", nrow(om_grid))

for (i in seq_len(nrow(om_grid))) {
  cat("Running OM", i, "/9 —", om_grid$om_id[i], "...")
  t0_i <- proc.time()
  mse_list[[i]] <- tryCatch(
    runMSE(OM = om_list[[i]], MPs = MPs_names, parallel = FALSE),
    error = function(e) {
      cat(" ERROR:", conditionMessage(e), "\n"); NULL
    }
  )
  dt <- round((proc.time() - t0_i)["elapsed"], 1)
  if (!is.null(mse_list[[i]])) cat(" OK (", dt, "s)\n") else cat(" FAILED\n")
}

n_ok <- sum(!sapply(mse_list, is.null))
cat("\n✅ MSE complete —", n_ok, "/", nrow(om_grid), "OMs successful\n")

if (n_ok == 0) stop("All OMs failed. Review errors.")

# Diagnostics for base OM (sigmaR=1.2, R0=0%)
idx_base <- which(om_grid$sigmaR == 1.2 & om_grid$R0_shift == 0)
mse_test  <- mse_list[[idx_base]]
if (!is.null(mse_test)) {
  cat("\n=== Base OM diagnostics (sigmaR=1.2, R0=0%) ===\n")
  for (j in seq_along(MPs_names)) {
    cat(MPs_names[j],
        "| year  1:", round(median(mse_test@SSB[, j,  1])),
        "| year 15:", round(median(mse_test@SSB[, j, 15])),
        "| year 30:", round(median(mse_test@SSB[, j, 30])), "t\n")
  }
}

# -------------------------------------------------------------
# STEP 10: PERFORMANCE METRICS
# -------------------------------------------------------------

metrics <- map_dfr(seq_len(nrow(om_grid)), function(i) {
  if (is.null(mse_list[[i]])) return(NULL)
  mse <- mse_list[[i]]
  map_dfr(seq_along(MPs_names), function(j) {
    B_rel <- mse@SSB[, j, ] / B0
    C     <- mse@Catch[, j, ]
    tibble(
      om_id            = om_grid$om_id[i],
      mp               = MPs_names[j],
      sigmaR           = om_grid$sigmaR[i],
      R0_shift         = om_grid$R0_shift[i],
      P_no_collapse    = mean(B_rel > 0.2),          # CCAMLR critical limit
      P_ccamlr_target  = mean(B_rel > 0.75),         # CCAMLR target
      mean_catch_kt    = mean(C) / 1e3,
      P_final_collapse = mean(B_rel[, n_proyear] < 0.2),
      P_final_target   = mean(B_rel[, n_proyear] > 0.75)
    )
  })
})

cat("\n=== Performance Metrics ===\n")
print(metrics |>
        select(om_id, mp, P_no_collapse, P_ccamlr_target,
               mean_catch_kt, P_final_collapse) |>
        mutate(across(where(is.numeric), ~round(., 3))))

# -------------------------------------------------------------
# STEP 11: FIGURES  —  all saved to fig/
# -------------------------------------------------------------

if (!dir.exists("fig")) dir.create("fig")

years_proj <- (yr_final + 1):(yr_final + n_proyear)
mp_labels  <- c("CC145" = "CC 145kt", "HCR_proportional" = "HCR Proportional")
mp_colors  <- c("CC 145kt" = "steelblue", "HCR Proportional" = "darkorange")

# ── Wormplot 3×3 ────────────────────────────────────────────

worm_data <- map_dfr(seq_len(nrow(om_grid)), function(i) {
  if (is.null(mse_list[[i]])) return(NULL)
  map_dfr(seq_along(MPs_names), function(j) {
    B_rel <- mse_list[[i]]@SSB[, j, ] / B0
    tibble(
      year     = years_proj,
      q05      = apply(B_rel, 2, quantile, 0.05),
      q25      = apply(B_rel, 2, quantile, 0.25),
      q50      = apply(B_rel, 2, quantile, 0.50),
      q75      = apply(B_rel, 2, quantile, 0.75),
      q95      = apply(B_rel, 2, quantile, 0.95),
      sigmaR   = om_grid$sigmaR[i],
      R0_shift = om_grid$R0_shift[i],
      MP_label = mp_labels[MPs_names[j]]
    )
  })
}) |>
  mutate(
    climate = factor(paste0("R0: ", R0_shift * 100, "%"),
                     levels = c("R0: 0%", "R0: -20%", "R0: -40%")),
    variab  = factor(paste0("σR=", sigmaR),
                     levels = c("σR=0.8", "σR=1.2", "σR=1.5"))
  )

p_worm <- ggplot(worm_data, aes(x = year)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = MP_label), alpha = 0.12) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = MP_label), alpha = 0.30) +
  geom_line(aes(y = q50, color = MP_label), linewidth = 0.9) +
  geom_hline(yintercept = 0.75, linetype = "dashed",
             color = "darkgreen", linewidth = 0.5) +
  geom_hline(yintercept = 0.20, linetype = "dashed",
             color = "red",       linewidth = 0.5) +
  facet_grid(variab ~ climate) +
  scale_color_manual(values = mp_colors, name = "HCR") +
  scale_fill_manual(values  = mp_colors, name = "HCR") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0, 2)) +
  geom_line(data        = hist_ssb,
            aes(x = year, y = ssb_rel),
            inherit.aes = FALSE,
            color       = "black",
            linewidth   = 0.7,
            linetype    = "solid") +
  geom_vline(xintercept = yr_final, linetype = "dotted",
             color = "gray40", linewidth = 0.4) +
  labs(x = "Year", y = "SSB / B0",
       title = "Relative Spawning Biomass Trajectories") +
  theme_bw() +
  theme(panel.grid.minor  = element_blank(),
        strip.background  = element_blank(),
        strip.text        = element_text(size = 9, face = "bold"),
        axis.text.x       = element_text(angle = 90, 
                                         size=10,
                                         hjust = 1),
        legend.position   = "bottom")

# ── Trade-off plot ───────────────────────────────────────────

clima_cols  <- c("0" = "steelblue", "-0.2" = "orange", "-0.4" = "red3")
clima_labs  <- c("0" = "No change", "-0.2" = "−20% R0", "-0.4" = "−40% R0")
mp_shapes   <- c("CC145" = 16, "HCR_proportional" = 17)
mp_labs2    <- c("CC145" = "CC 145kt", "HCR_proportional" = "HCR Proportional")

p_tradeoff <- ggplot(metrics,
                     aes(x = mean_catch_kt,
                         y = P_ccamlr_target,
                         color = as.character(R0_shift),
                         shape = mp)) +
  geom_point(size = 4.5, alpha = 0.9) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "darkgreen", alpha = 0.5) +
  scale_color_manual(values = clima_cols, labels = clima_labs, name = "Climate scenario") +
  scale_shape_manual(values = mp_shapes,  labels = mp_labs2,   name = "HCR") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Mean catch (thousand t)", y = "P(SSB > 0.75 B0)",
       title = "Trade-off: Yield vs. CCAMLR Target") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

# ── Performance metrics barplot ──────────────────────────────

metrics_long <- metrics |>
  group_by(mp) |>
  summarise(
    `P(B > 0.75B0)` = mean(P_ccamlr_target),
    `P(B > 0.20B0)` = mean(P_no_collapse),
    `P(col. t30)`   = mean(P_final_collapse),
    `Mean catch`    = mean(mean_catch_kt) / max(mean(mean_catch_kt))  # normalized 0-1
  ) |>
  pivot_longer(-mp, names_to = "metric", values_to = "value") |>
  mutate(MP_label = mp_labs2[mp])

p_barras <- ggplot(metrics_long, aes(x = metric, y = value, fill = MP_label)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.65) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = mp_colors, name = "HCR") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = NULL, y = "Probability / relative value",
       title = "Average performance metrics by HCR") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x      = element_text(size = 8),
        legend.position  = "bottom")

# ── CDF of SSB/B0 at final year ──────────────────────────────

cdf_data <- map_dfr(seq_len(nrow(om_grid)), function(i) {
  if (is.null(mse_list[[i]])) return(NULL)
  map_dfr(seq_along(MPs_names), function(j) {
    B_rel_final <- mse_list[[i]]@SSB[, j, n_proyear] / B0
    tibble(B_rel = B_rel_final,
           MP_label = mp_labels[MPs_names[j]])
  })
})

p_cdf <- ggplot(cdf_data, aes(x = B_rel, color = MP_label)) +
  stat_ecdf(linewidth = 0.9) +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "darkgreen", linewidth = 0.5) +
  geom_vline(xintercept = 0.20, linetype = "dashed", color = "red",       linewidth = 0.5) +
  scale_color_manual(values = mp_colors, name = "HCR") +
  scale_x_continuous(limits = c(0, 2.5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "SSB / B0 (year 30)", y = "Cumulative probability",
       title = "CDF of relative biomass — final year") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# ── Catch boxplot ────────────────────────────────────────────

catch_data <- map_dfr(seq_len(nrow(om_grid)), function(i) {
  if (is.null(mse_list[[i]])) return(NULL)
  map_dfr(seq_along(MPs_names), function(j) {
    tibble(
      catch_kt = as.vector(mse_list[[i]]@Catch[, j, ]) / 1e3,
      MP_label = mp_labels[MPs_names[j]],
      climate  = clima_labs[as.character(om_grid$R0_shift[i])]
    )
  })
}) |> mutate(climate = factor(climate, levels = c("No change", "−20% R0", "−40% R0")))

p_catch <- ggplot(catch_data, aes(x = climate, y = catch_kt, fill = MP_label)) +
  geom_boxplot(alpha = 0.8, outliers = FALSE) +
  scale_fill_manual(values = mp_colors, name = "HCR") +
  labs(x = NULL, y = "Annual catch (thousand t)",
       title = "Catch distribution by climate scenario") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# ── Recruitment wormplot ────────────────────────────────────

# Extract recruitment trajectories (age 0 numbers from N_ageArray)
# N_ageArray dims: [nsim, n_mp, n_age, n_yr_proj]
# Age 0 = first column of age dimension
# ── Summary table ────────────────────────────────────────────

tabla_df <- metrics |>
  mutate(
    MP            = mp_labels[mp],
    `R0 shift`    = paste0(R0_shift * 100, "%"),
    `P(B>0.2B0)`  = scales::percent(P_no_collapse,    accuracy = 0.1),
    `P(B>0.75B0)` = scales::percent(P_ccamlr_target,  accuracy = 0.1),
    `Catch(kt)`   = round(mean_catch_kt, 0),
    `P(col.t30)`  = scales::percent(P_final_collapse, accuracy = 0.1)
  ) |>
  select(sigmaR, `R0 shift`, MP, `P(B>0.2B0)`,
         `P(B>0.75B0)`, `Catch(kt)`, `P(col.t30)`)

p_tabla <- gridExtra::tableGrob(
  tabla_df, rows = NULL,
  theme = gridExtra::ttheme_minimal(base_size = 7.5)
)

# ── Combined main figure ─────────────────────────────────────

fig_main <- p_worm / (p_tradeoff | wrap_elements(p_tabla)) +
  plot_layout(heights = c(2.2, 1)) +
  plot_annotation(
    title = "Antarctic Krill MSE — HCR Comparison",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# ── Additional MSE figures ───────────────────────────────────

fig_extra <- (p_barras | p_cdf) / p_catch +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Antarctic Krill MSE — Performance Metrics",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# ── Save to fig/ ─────────────────────────────────────────────

ggsave("fig/MSE_krill_openMSE.png",  fig_main,  width = 18, height = 15, dpi = 150)
ggsave("fig/MSE_krill_metrics.png",  fig_extra, width = 16, height = 12, dpi = 150)

# Individual figures for RMD report
ggsave("fig/fig_wormplot.png",   p_worm,     width = 14, height = 10, dpi = 150)
#ggsave("fig/fig_recruitment.png", p_recruit,  width = 14, height = 10, dpi = 150)
ggsave("fig/fig_tradeoff.png",   p_tradeoff, width = 10, height =  7, dpi = 150)
ggsave("fig/fig_barplot.png",    p_barras,   width =  9, height =  6, dpi = 150)
ggsave("fig/fig_cdf.png",        p_cdf,      width =  9, height =  6, dpi = 150)
ggsave("fig/fig_catch.png",      p_catch,    width =  9, height =  6, dpi = 150)

cat("\n✅ fig/MSE_krill_openMSE.png\n")
cat("✅ fig/MSE_krill_metrics.png\n")
cat("✅ fig/fig_wormplot.png\n")
#cat("✅ fig/fig_recruitment.png\n")
cat("✅ fig/fig_tradeoff.png\n")
cat("✅ fig/fig_barplot.png\n")
cat("✅ fig/fig_cdf.png\n")
cat("✅ fig/fig_catch.png\n")
cat("✅ Script v7.0 completed successfully\n")

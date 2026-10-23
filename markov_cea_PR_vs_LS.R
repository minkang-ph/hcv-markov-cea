# =================================================================
# Markov cohort model - cost-effectiveness analysis (CEA)
# Chronic hepatitis C: PR (pegylated interferon + ribavirin)
#                 vs   LS (ledipasvir + sofosbuvir, the new regimen)
#
# Author: Min Kang
#
# Note: Model structure and input parameters are adapted from
# coursework materials, which draw on published studies.
# =================================================================


# -----------------------------------------------------------------
# STEP 1. INPUT PARAMETERS
# -----------------------------------------------------------------
# Every input is typed ONCE here; later steps read it from these lists.

params <- list(
  # --- probability of SVR (treatment success) ---
  p_svr_PR = 0.395,     # PR
  p_svr_LS = 0.89,      # LS

  # --- annual transition probabilities ---
  p_mod_cirr = 0.043,   # Moderate  -> Cirrhosis
  p_cirr_lf  = 0.039,   # Cirrhosis -> Liver failure
  p_lf_death = 0.135,   # Liver failure -> Death (disease-related)

  # --- cost of ONE cycle spent in each state ---
  cost_svr   = 0,       # assumed 0 (cured)
  cost_mod   = 2000,
  cost_cirr  = 3700,
  cost_lf    = 6300,
  cost_death = 0,       # cost of the transition Liver failure -> Death

  # --- quality-of-life weights (utilities) for one cycle ---
  u_svr  = 0.71,
  u_mod  = 0.66,
  u_cirr = 0.55,
  u_lf   = 0.45,

  # --- model settings ---
  start_age   = 40,
  n_cycles    = 70,
  cohort_size = 1000,
  disc_cost   = 0.035,  # discount rate for costs
  disc_qaly   = 0.035   # discount rate for QALYs and life years
)

# Standard errors (SE) of the uncertain parameters.
# Sensitivity analyses (mean +/- 1.96 x SE).
params_se <- list(
  p_svr_PR   = 0.02,
  p_svr_LS   = 0.03,
  p_mod_cirr = 0.01,
  p_cirr_lf  = 0.01,
  p_lf_death = 0.02,
  cost_mod   = 230,
  cost_cirr  = 1100,
  cost_lf    = 1900,
  u_svr      = 0.052,
  u_mod      = 0.054,
  u_cirr     = 0.054,
  u_lf       = 0.035
)

# Age-specific natural mortality rates
#   0-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+
age_breaks <- c(0,     30,    40,    50,    60,    70,    80)
mort_rates <- c(0.002, 0.002, 0.004, 0.012, 0.027, 0.066, 0.189)

# Cost of the full treatment course per patient, paid once at cycle 0 (the start)
cost_tx_PR <- 6400
cost_tx_LS <- 44000


# -----------------------------------------------------------------
# STEP 2. HELPER FUNCTION: the transition matrix for ONE cycle
# -----------------------------------------------------------------
# A transition matrix says: "of the people in state X (row), what share
# move to state Y (column) during one cycle?"  Each row must add up to 1.
# mu = natural death risk at the current age.

state_names <- c("SVR", "Moderate", "Cirrhosis", "LiverFailure", "Death")

build_transition_matrix <- function(mu, params) {

  P <- matrix(0, nrow = 5, ncol = 5,
              dimnames = list(from = state_names, to = state_names))

  # SVR: cured, can only die of natural causes
  P["SVR", "Death"] <- mu
  P["SVR", "SVR"]   <- 1 - mu

  # Moderate: may progress to cirrhosis, or die naturally
  P["Moderate", "Cirrhosis"] <- params$p_mod_cirr
  P["Moderate", "Death"]     <- mu
  P["Moderate", "Moderate"]  <- 1 - params$p_mod_cirr - mu

  # Cirrhosis: may progress to liver failure, or die naturally
  P["Cirrhosis", "LiverFailure"] <- params$p_cirr_lf
  P["Cirrhosis", "Death"]        <- mu
  P["Cirrhosis", "Cirrhosis"]    <- 1 - params$p_cirr_lf - mu

  # Liver failure: dies of the disease OR naturally
  P["LiverFailure", "Death"]        <- mu + params$p_lf_death
  P["LiverFailure", "LiverFailure"] <- 1 - mu - params$p_lf_death

  # Death is the absorbing state: once dead, always dead
  P["Death", "Death"] <- 1

  # Safety check: every row must sum to 1, otherwise stop with an error
  stopifnot(all(abs(rowSums(P) - 1) < 1e-10))

  return(P)
}


# -----------------------------------------------------------------
# STEP 3. THE MARKOV MODEL FUNCTION
# -----------------------------------------------------------------
run_markov <- function(p_svr, cost_tx, params) {

  n <- params$n_cycles

  # The "Markov trace": rows = cycles, columns = states.
  trace <- matrix(0, nrow = n, ncol = 5,
                  dimnames = list(NULL, state_names))

  # Cycle 1 = initial split of the cohort by treatment success.
  trace[1, "SVR"]      <- params$cohort_size * p_svr
  trace[1, "Moderate"] <- params$cohort_size * (1 - p_svr)

  # Cycles 2..n: multiply last cycle's distribution by the transition matrix
  for (t in 2:n) {
    age <- params$start_age + t                       # age in this cycle
    mu  <- mort_rates[findInterval(age, age_breaks)]  # natural death risk
    P   <- build_transition_matrix(mu, params)
    trace[t, ] <- trace[t - 1, ] %*% P                # matrix multiplication
  }

  # ---- outcomes per cycle ----
  cycle  <- 1:n
  age    <- params$start_age + cycle
  alive  <- rowSums(trace[, 1:4])                     # everyone not in Death

  # QALYs = people in each state x that state's quality-of-life weight
  u_vec <- c(params$u_svr, params$u_mod, params$u_cirr, params$u_lf, 0)
  qaly  <- as.vector(trace %*% u_vec)

  # Costs = people in each state x that state's cost per cycle
  c_vec <- c(params$cost_svr, params$cost_mod, params$cost_cirr,
             params$cost_lf, 0)

  # New liver-failure deaths this cycle = people in liver failure LAST cycle
  # x death probability (c(0, head(..., -1)) shifts the column down one cycle)
  new_lf_deaths <- c(0, head(trace[, "LiverFailure"], -1)) * params$p_lf_death

  # Cost of the transition to death is 0 in the base case; it can be varied
  # in a sensitivity analysis
  cost_cyc <- as.vector(trace %*% c_vec) + new_lf_deaths * params$cost_death

  # ---- discounting: 1 / (1 + rate)^year ----
  df_cost <- 1 / (1 + params$disc_cost) ^ cycle
  df_out  <- 1 / (1 + params$disc_qaly) ^ cycle   
    # for life-years and QALYs

  trace_df <- data.frame(
    cycle, age,
    round(trace, 2),  # n in each state, rounded for display only
    total_death = round(params$cohort_size - alive, 2),
    ly_undisc   = alive,
    qaly_undisc = qaly,
    cost_undisc = cost_cyc,
    df_cost, df_out,
    ly_disc     = alive * df_out,
    qaly_disc   = qaly * df_out,
    cost_disc   = cost_cyc * df_cost
  )

  # ---- lifetime totals for the whole cohort ----
  # Treatment cost is paid once at the start = no discount
  total_cost <- cost_tx * params$cohort_size + sum(trace_df$cost_disc)

  summary <- data.frame(
    ly_undisc   = sum(trace_df$ly_undisc),
    ly_disc     = sum(trace_df$ly_disc),
    qaly_undisc = sum(trace_df$qaly_undisc),
    qaly_disc   = sum(trace_df$qaly_disc),
    cost_disc   = total_cost
  )

  # Wrap the two tables in a list so the function can return both
  return(list(trace = trace_df, summary = summary))
}


# -----------------------------------------------------------------
# STEP 4. RUN THE MODEL FOR BOTH STRATEGIES & CHECK
# -----------------------------------------------------------------
PR <- run_markov(params$p_svr_PR, cost_tx_PR, params)
LS <- run_markov(params$p_svr_LS, cost_tx_LS, params)

cat("\n--- Sanity checks ---\n")

# Starting split check
print(PR$trace[1, c("SVR", "Moderate")])   # expect 395 SVR / 605 Moderate
print(LS$trace[1, c("SVR", "Moderate")])   # expect 890 SVR / 110 Moderate

# Every cycle, all 5 states should add up to the 1000 patients (expect TRUE)
# Tolerance 0.03 > max rounding error of 5 states x 0.005 = 0.025
print(all(abs(rowSums(PR$trace[, c("SVR","Moderate","Cirrhosis","LiverFailure","Death")]) - 1000) < 0.03))
print(all(abs(rowSums(LS$trace[, c("SVR","Moderate","Cirrhosis","LiverFailure","Death")]) - 1000) < 0.03))

# By the last cycle nearly everyone should have died (expect about 1000)
print(tail(PR$trace$Death, 1))
print(tail(LS$trace$Death, 1))

# Totals for the whole cohort
cat("\n--- Cohort totals (", params$cohort_size, "patients ) ---\n")
cat("PR:\n")
print(round(PR$summary, 1))
cat("LS:\n")
print(round(LS$summary, 1))


# -----------------------------------------------------------------
# STEP 5. COST-EFFECTIVENESS RESULTS (per patient)
# -----------------------------------------------------------------
n_pat <- params$cohort_size

results <- data.frame(
  Strategy = c("PR", "LS"),
  Cost     = c(PR$summary$cost_disc,   LS$summary$cost_disc)   / n_pat,
  LYs      = c(PR$summary$ly_disc,     LS$summary$ly_disc)     / n_pat,
  QALYs    = c(PR$summary$qaly_disc,   LS$summary$qaly_disc)   / n_pat
)

cat("\n--- Discounted results per patient ---\n")
print(round(results[, -1], 3))   # all columns except the first (Strategy name)

# Incremental analysis: LS minus PR
d_cost <- results$Cost[2]  - results$Cost[1]
d_qaly <- results$QALYs[2] - results$QALYs[1]
d_ly   <- results$LYs[2]   - results$LYs[1]

cat("\nIncremental cost  :", round(d_cost, 0), "\n")
cat("Incremental QALYs :", round(d_qaly, 3), "\n")
cat("Incremental LYs   :", round(d_ly, 3), "\n")

# ICER = extra cost / extra QALYs
if (d_qaly > 0 && d_cost < 0) {
  cat("LS DOMINATES PR\n")
} else if (d_qaly > 0) {
  cat("ICER:", round(d_cost / d_qaly, 0), "per QALY gained\n")
} else {
  cat("LS gives no QALY gain over PR - check inputs.\n")
}


# -----------------------------------------------------------------
# STEP 6. DRAW A PLOT
# -----------------------------------------------------------------
cols <- c("darkgreen", "orange", "red", "purple", "black")

par(mfrow = c(1, 2))          # two plots side by side

# Left plot: PR
matplot(PR$trace$age, PR$trace[, state_names],
        type = "l", lty = 1, lwd = 2, col = cols,
        xlab = "Age", ylab = "Patients", main = "Markov trace: PR")
legend("right", legend = state_names, col = cols, lty = 1, lwd = 2,
       bty = "n", cex = 0.7)

# Right plot: LS
matplot(LS$trace$age, LS$trace[, state_names],
        type = "l", lty = 1, lwd = 2, col = cols,
        xlab = "Age", ylab = "Patients", main = "Markov trace: LS")
legend("right", legend = state_names, col = cols, lty = 1, lwd = 2,
       bty = "n", cex = 0.7)

par(mfrow = c(1, 1))          # settings back to one plot at a time


# -----------------------------------------------------------------
# STEP 7. DETERMINISTIC (ONE-WAY) SENSITIVITY ANALYSIS
# -----------------------------------------------------------------
# Identifies which parameters have the most impact on the NMB (and ICER) of LS vs PR.
# Each parameter with a standard error is varied one at a time across its 95% range ; the discount rate is varied as a scenario.
# Also includes a threshold analysis for the maximum price of LS.

# --- 7a. Willingness-to-pay (WTP) threshold ---
gdp_pc   <- 77350.76        # GDP per capita 2023 Egypt
wtp_base <- 0.39 * gdp_pc   # 39% of GDP per capita (Riviere) = about 30,167

# --- 7b. One function that runs the whole CEA for a given set of inputs ---
# inp   = list of all inputs, same structure with params
# tx_LS = price of LS (default = base-case price); varied in the threshold analysis
# NMB = incremental QALYs x WTP per QALY - incremental cost (whole cohort)

run_cea <- function(inp, wtp, tx_LS = cost_tx_LS) {
  PR_res <- run_markov(inp$p_svr_PR, cost_tx_PR, inp) # fixed price for PR
  LS_res <- run_markov(inp$p_svr_LS, tx_LS, inp) # varied price for LS
  # run_markov <- function(p_svr, cost_tx, params)
  
  d_qaly <- LS_res$summary$qaly_disc - PR_res$summary$qaly_disc
  d_cost <- LS_res$summary$cost_disc - PR_res$summary$cost_disc
  nmb    <- d_qaly * wtp - d_cost

  return(c(d_qaly = d_qaly, d_cost = d_cost, icer = d_cost / d_qaly, nmb = nmb))
}

# Helper: change ONE input by name. The discount rate changes both rates.
# inp[[name]] picks the item whose name is stored in the variable 'name'.
set_param <- function(inp, name, value) {
  if (name == "discount_rate") {
    inp$disc_cost <- value
    inp$disc_qaly <- value 
    # Cost and QALY discount rates are varied together (same value), by design choice
  } else {
    inp[[name]] <- value  
    # in inp, find the item called <name> (e.g. "u_svr") and set it to <value>
  }
  return(inp)
}

base_inputs <- params   # all inputs, including the SVR rates, are in params

base_res <- run_cea(base_inputs, wtp_base)
cat("\nBase case: ICER =", round(base_res[["icer"]]),
    "| NMB =", round(base_res[["nmb"]]), "\n")

# --- 7c. Table of parameters with their bounds (mean +/- 1.96 x SE) ---
# Note: at its low bound the SVR utility (0.608) falls below the Moderate utility (0.66) - see README limitations.

dsa_table <- data.frame(
  name = names(params_se),  
  # the 12 item names of params_se, as a text (character) vector
  mean = unlist(params[names(params_se)], use.names = FALSE), 
  # base values from params, dropping the names
  se   = unlist(params_se, use.names = FALSE), 
  stringsAsFactors = FALSE # use names instead of positions (factors)
)
dsa_table$lower <- dsa_table$mean - 1.96 * dsa_table$se
dsa_table$upper <- dsa_table$mean + 1.96 * dsa_table$se

# Discount rate: no sampling error (design choice), so a scenario range of 2% to 6% instead of SE
dsa_table <- rbind(dsa_table, # Append discount-rate row underneath (rbind requires the same column names)
                   data.frame(name = "discount_rate", mean = params$disc_cost,
                              se = NA, lower = 0.02, upper = 0.06))

# Display labels for the tornado diagram
param_labels <- c(
  p_svr_PR      = "Probability - SVR in PR",
  p_svr_LS      = "Probability - SVR in LS",
  p_mod_cirr    = "Probability - From M to C",
  p_cirr_lf     = "Probability - From C to LF",
  p_lf_death    = "Probability - From LF to Death",
  cost_mod      = "M costs",
  cost_cirr     = "C costs",
  cost_lf       = "LF costs",
  u_svr         = "SVR QoL",
  u_mod         = "M QoL",
  u_cirr        = "C QoL",
  u_lf          = "LF QoL",
  discount_rate = "Discount Rate"
)
dsa_table$label <- unname(param_labels[dsa_table$name])
 # add a new column 'label' in dsa_table 
 # find the labels matching dsa_table$name, remove the tags, keep only the text

# --- 7d. Run the model at the low and high value of each parameter ---
dsa_results <- data.frame()   # start with an empty table

for (i in 1:nrow(dsa_table)) {
  low_res  <- run_cea(set_param(base_inputs, dsa_table$name[i], dsa_table$lower[i]), wtp_base)
  # For parameter i of dsa_table, replace its value in base_inputs with dsa_table$lower (in a copy, so base_inputs itself is unchanged), then run the CEA.
  high_res <- run_cea(set_param(base_inputs, dsa_table$name[i], dsa_table$upper[i]), wtp_base)
  #run_cea <- function(inp, wtp, tx_LS = cost_tx_LS) 
  # no SE (statistical uncertainty) for tx prices; threshold analysis

  # add one row for this parameter
  dsa_results <- rbind(dsa_results, data.frame(
    parameter  = dsa_table$label[i],
    dqaly_low  = low_res[["d_qaly"]], dqaly_high = high_res[["d_qaly"]],
    dcost_low  = low_res[["d_cost"]], dcost_high = high_res[["d_cost"]],
    nmb_low    = low_res[["nmb"]],    nmb_high   = high_res[["nmb"]],
    icer_low   = low_res[["icer"]],   icer_high  = high_res[["icer"]]
  ))
}

# Size of the NMB swing, then sort biggest first
dsa_results$nmb_range <- abs(dsa_results$nmb_high - dsa_results$nmb_low)
  # add new column nmb_range to the dataframe
dsa_results <- dsa_results[order(dsa_results$nmb_range, decreasing = TRUE), ]

cat("\n--- One-way sensitivity analysis (cohort of", params$cohort_size, ") ---\n")
dsa_print <- dsa_results
dsa_print[, -1] <- round(dsa_print[, -1], 0)   
  # make a copy of dsa_results to round up all column values, except the name
print(dsa_print, row.names = FALSE)  # hide automatic row numbers

# --- 7e. Tornado diagram (orange = low input value, green = high input value) ---
plot_tornado <- function(res, base_nmb) {
  res  <- res[order(res$nmb_range), ]   
  # smallest from y axis=1, so the biggest ends up on top
  n    <- nrow(res)
  lo   <- res$nmb_low  / 1e6            # in millions, easier to read
  hi   <- res$nmb_high / 1e6
  base <- base_nmb / 1e6

  old_par <- par(mar = c(5, 13, 4, 2))  
  # wide left margin for the long labels on the left
  # c(bottom, left, top, right)
  
  plot(NA, xlim = range(c(lo, hi, base, 0)), ylim = c(0.5, n + 0.5),
       yaxt = "n", ylab = "", # no default y axis & drop its title
       xlab = "Net monetary benefit (millions, cohort of 1,000)",
       main = "Tornado diagram: NMB of LS vs PR")

  for (i in 1:n) {
    # bar from base NMB to the NMB with parameter i at its LOW values;
    # 0.8 thick; goes left if that NMB < base NMB, right if it is higher
    rect(base, i - 0.4, lo[i], i + 0.4, col = "darkorange")
    # same, for HIGH values
    rect(base, i - 0.4, hi[i], i + 0.4, col = "darkgreen")
  }
  axis(2, at = 1:n, labels = res$parameter, las = 1, cex.axis = 0.8)
    #draw axis on the left            # horizontal text
                                              # 80% character expansion
  abline(v = base, lwd = 2)             
    # vertical line at base-case NMB, twice the width
  abline(v = 0, lty = 2, col = "red")   
   # NMB = 0, dashed line
   # left of this, LS is NOT cost-effective
  legend("bottomright", legend = c("Low input value", "High input value"),
         fill = c("darkorange", "darkgreen"), bty = "n", cex = 0.8)
                                        # box type=none
  par(old_par)              # restore the original margins
}

plot_tornado(dsa_results, base_res[["nmb"]])

# --- 7f. # Threshold (maximum acceptable) price analysis ---

# NMB is linear in the LS price P:  NMB(LS price of P) = NMB(LS price of 0) - cohort size * LS price of P
# Setting NMB(P*) = 0 gives the threshold price P* = NMB(0) / 1,000.
# LS is cost-effective at the given WTP for any price P <= P*.

wtp_values <- c(Low = 0.20, Base = 0.39, High = 0.49) * gdp_pc
  # source of range: Riviere
  # assigned name & values in wtp_values (named vector)

cat("\n--- Maximum LS price per patient ---\n")
for (i in 1:3) { # for three different wtp_values (low, base, high)
  nmb_at_zero_price <- run_cea(base_inputs, wtp_values[[i]], tx_LS = 0)[["nmb"]] 
  # run CEA with LS price = 0, keep only the nmb value, save it as nmb_at_zero_price
  max_price <- nmb_at_zero_price / params$cohort_size
  cat(names(wtp_values)[i], "WTP =", round(wtp_values[[i]]),
      "-> max LS price =", round(max_price), "\n")
} # Maximum cost effective price for each wtp thresholds


# -----------------------------------------------------------------
# STEP 8. PROBABILISTIC SENSITIVITY ANALYSIS (Monte Carlo simulation)
# -----------------------------------------------------------------
# Only inputs with a standard error (sampling uncertainty) are drawn from distributions: SVR rates, transition probabilities and utilities (beta), and state costs (gamma).
# Kept fixed: discount rate (methodological choice), treatment prices (known values), natural mortality (life-table values).


# --- 8a. Which distribution for which parameter ---
# Probabilities and utilities lie between 0 and 1 -> beta distribution
# Costs are positive and right-skewed             -> gamma distribution

beta_names  <- c("p_svr_PR", "p_svr_LS", "p_mod_cirr", "p_cirr_lf", "p_lf_death",
                 "u_svr", "u_mod", "u_cirr", "u_lf")
gamma_names <- c("cost_mod", "cost_cirr", "cost_lf")

# --- 8b. Function: draw ONE random set of inputs ---
# Beta (values between 0 and 1): from mean m and SE s 
#   common = alpha + beta = m(1-m)/s^2 - 1
#   alpha = m * common;  beta = (1-m) * common     (R: shape1 = alpha, shape2 = beta)
# Gamma (positive, right-skewed values): 
# shape = m^2 / s^2;  scale = s^2 / m

draw_params <- function(inp, se) {
  for (nm in beta_names) {
    m      <- inp[[nm]]                 # mean
    s      <- se[[nm]]                  # standard error
    common <- m * (1 - m) / s^2 - 1
    inp[[nm]] <- rbeta(1, shape1 = m * common, shape2 = (1 - m) * common) # draw 1 random number from beta distribution (shape 1, shape 2), 
    # replace the value in the copy of inp (original params stay)
  }
  for (nm in gamma_names) {
    m <- inp[[nm]]
    s <- se[[nm]]
    inp[[nm]] <- rgamma(1, shape = m^2 / s^2, scale = s^2 / m)
  }
  return(inp)
}

# --- 8c. Monte Carlo simulation ---
n_sim <- 1000
set.seed(123)    # makes the random draws reproducible

# empty results table with one row per simulation
psa <- data.frame(d_qaly = numeric(n_sim), d_cost = numeric(n_sim))
 # numeric vector of length n_sim filled entirely with zeros

for (sim in 1:n_sim) {
  draw <- draw_params(base_inputs, params_se)   # one random set of inputs
  res  <- run_cea(draw, wtp_base)               # run PR vs LS with that set
  psa$d_qaly[sim] <- res[["d_qaly"]]            # store incremental QALYs # could have also said psa[sim, "d_qaly"]
  psa$d_cost[sim] <- res[["d_cost"]]            # store incremental cost
}

cat("\n--- PSA results (", n_sim, "simulations, cohort of", params$cohort_size, ") ---\n")
cat("Mean incremental QALYs :", round(mean(psa$d_qaly), 1), "\n")
cat("Mean incremental cost  :", round(mean(psa$d_cost)), "\n")
cat("ICER of the means      :", round(mean(psa$d_cost) / mean(psa$d_qaly)), "\n")

# Probability of being cost-effective at the three WTP thresholds
for (i in 1:3) {
  nb <- psa$d_qaly * wtp_values[[i]] - psa$d_cost   # net benefit in each simulation
  cat(names(wtp_values)[i], "WTP =", round(wtp_values[[i]]),
      "-> probability cost-effective =", round(mean(nb > 0), 3), "\n")
}
  # count nb>0 TRUE values in 1,000 rows and round the number

# --- 8d. Cost-effectiveness plane ---
# x axis: QALY ; y axis: costs (per thousand for presentation)

wtp_cols <- c("red", "darkorange", "darkgreen")   # 20%, 39%, 49% of GDP

old_par <- par(mar = c(5, 3, 4, 2))
               
plot(psa$d_qaly, psa$d_cost / 1000,
 pch = 16, cex = 0.5, col = "navy",
 xlab = "Incremental effectiveness (QALYs, cohort of 1,000)",
 ylab = "Incremental cost (thousands)",
 main = "Cost-effectiveness plane: LS vs PR")
               
 abline(h = 0, v = 0, col = "grey")
               
for (i in 1:3) {
 abline(a = 0, b = wtp_values[[i]] / 1000,
 col = wtp_cols[i], lwd = 2)
}
               
legend("bottomright", legend = c("20% GDP", "39% GDP", "49% GDP"),
col = wtp_cols, lwd = 2, bty = "n")
               
               
par(old_par) # Restore original graphical settings

# --- 8e. Cost-effectiveness acceptability curve (CEAC) ---
wtp_grid <- seq(0, 50000, by = 1000) # 0, 1000, 2000, 3000 ...
ceac <- data.frame(wtp = wtp_grid, prob_ce = NA)
 # make an empty data frame, assign two columns wtp and prob_ce

for (k in 1:length(wtp_grid)) {
  nb <- psa$d_qaly * wtp_grid[k] - psa$d_cost    # net benefit at this WTP, per simulation
  ceac$prob_ce[k] <- mean(nb > 0)  # share of simulations with net benefit > 0
}

old_par <- par(mar = c(5, 3, 4, 2))

plot(ceac$wtp, ceac$prob_ce, type = "l", lwd = 2, col = "darkorange",
     ylim = c(0, 1), # y axis limits from 0 to 1
     xlab = "Willingness-to-pay threshold (per QALY)",
     ylab = "Probability cost-effective",
     main = "Cost-effectiveness acceptability curve: LS vs PR")
abline(v = wtp_values, lty = 3, col = wtp_cols)  # the three WTP thresholds

par(old_par)

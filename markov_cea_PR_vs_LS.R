#=================================================================
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


params <- list(
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
  disc_cost   = 0.035,   # discount rate for costs
  disc_qaly   = 0.035  # discount rate for QALYs and life years
)

# Age specific mortality rates
#   0-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+
age_breaks <- c(0,     30,    40,    50,    60,    70,    80)
mort_rates <- c(0.002, 0.002, 0.004, 0.012, 0.027, 0.066, 0.189)

# Treatment-specific inputs
p_svr_PR <- 0.395    # probability of SVR with PR
p_svr_LS <- 0.89     # probability of SVR with LS

# Cost of the full treatment course per patient, paid once at cycle 0 (the start).
cost_tx_PR <- 6400
cost_tx_LS <- 44000


# -----------------------------------------------------------------
# STEP 2. HELPER FUNCTION: the transition matrix for ONE cycle
# -----------------------------------------------------------------
# A transition matrix says: "of the people in state X (row), what share move to state Y (column) during one cycle?"  Each row must add up to 1.
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
    age <- params$start_age + t   # age in this cycle
    mu  <- mort_rates[findInterval(age, age_breaks)]  # natural death risk
    P   <- build_transition_matrix(mu, params)
    trace[t, ] <- trace[t - 1, ] %*% P   # matrix multiplication
  }
  
  # ---- outcomes per cycle ----
  cycle  <- 1:n
  age    <- params$start_age + cycle
  alive  <- rowSums(trace[, 1:4])  # everyone not in Death
  ly     <- alive          # life-years this cycle
  
  # QALYs = people in each state x that state's quality-of-life weight
  u_vec    <- c(params$u_svr, params$u_mod, params$u_cirr, params$u_lf, 0)
  qaly     <- as.vector(trace %*% u_vec)
  
  # Costs = people in each state x that state's cost per cycle
  c_vec    <- c(params$cost_svr, params$cost_mod, params$cost_cirr,
                params$cost_lf, 0)
  new_lf_deaths <- c(0, head(trace[, "LiverFailure"], -1)) * params$p_lf_death
  # New liver-failure deaths this cycle = people in liver failure LAST cycle x death probability
  # (c(0, head(..., -1)) shifts the column down by one cycle)

    cost_cyc <- as.vector(trace %*% c_vec) + new_lf_deaths * params$cost_death
    #sensitivity analysis of cost of death !=0
  
  # ---- discounting: 1 / (1 + rate)^year ----
  # Life-years and QALYs discount rate: df_out.
  df_cost <- 1 / (1 + params$disc_cost) ^ cycle
  df_out  <- 1 / (1 + params$disc_qaly) ^ cycle

  
  trace_df <- data.frame(
    cycle, age, round(trace, 2), # n(each state) rounded to 2nd decimal for display
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
  
  return(list(trace = trace_df, summary = summary))
  # Wrap the two tables in a list so the function can return both
}

# -----------------------------------------------------------------
# STEP 4. RUN THE MODEL FOR BOTH STRATEGIES & CHECK
# -----------------------------------------------------------------
PR <- run_markov(p_svr_PR, cost_tx_PR, params)
LS <- run_markov(p_svr_LS, cost_tx_LS, params)

# First check: PR/LS totals per 1000 patients. 
cat("\n--- PR/LS totals for the cohort of", params$cohort_size, "---\n")

# Starting Split check
PR$trace[1, c("SVR", "Moderate")] # 395 SVR / 605 Moderate
LS$trace[1, c("SVR", "Moderate")] # 890 SVR / 110 Moderate

# Every cycle, all 5 states should add up to the 1000 patients
# Tolerance 0.03 > max rounding error of 5 states x 0.005 = 0.025
all(abs(rowSums(PR$trace[, c("SVR","Moderate","Cirrhosis","LiverFailure","Death")]) - 1000) < 0.03)
all(abs(rowSums(LS$trace[, c("SVR","Moderate","Cirrhosis","LiverFailure","Death")]) - 1000) < 0.03)

# By the last cycle nearly everyone should have died
tail(PR$trace$Death, 1)
tail(LS$trace$Death, 1)

print(round(PR$summary, 1))
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
print(round(results[, -1], 3)) # except for the first column (Strategy name)

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
legend("right", legend = state_names, col = cols, lty = 1, lwd = 2, bty = "n", cex = 0.7)

# Right plot: LS
matplot(LS$trace$age, LS$trace[, state_names],
        type = "l", lty = 1, lwd = 2, col = cols,
        xlab = "Age", ylab = "Patients", main = "Markov trace: LS")
legend("right", legend = state_names, col = cols, lty = 1, lwd = 2, bty = "n", cex = 0.7)

par(mfrow = c(1, 1))          # settings back to one plot at a time
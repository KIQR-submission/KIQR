# Title: Main Simulation Execution Script
# Description: This script configures and runs the KIQR simulation study.
#              It sources the core package functions.

# --- 1. SETUP ----
# Load all libraries and source all necessary functions
cat("Loading required packages and functions...\n")
library(MASS)
library(mvtnorm)
library(progress)
library(parallel)
library(quantreg)
library(rqPen)
library(dplyr)

# Source the refactored package functions
source("R/loss_penalty_functions.R")
source("R/tuning_criteria.R")
source("R/solvers.R")
source("R/kiqr_api.R")
source("R/utils.R")

# HELPER FUNCTION (only for the baseline rq.pen model)
run_pen_qr_replicate <- function(i, X_input, Y_input, tau, penalty_factor = NULL) {
  if (is.null(penalty_factor)) {
    penalty_factor <- rep(1, NCOL(X_input))
  }
  fit <- rqPen::rq.pen(
    x = X_input,
    y = as.vector(Y_input[i, ]),
    tau = tau,
    penalty.factor = penalty_factor,
    penalty = "SCAD"
  )
  best_lambda <- select_lambda_hbic(fit)
  as.vector(stats::coef(fit, lambda = best_lambda))
}


# --- 2. CONFIGURATION ----
cat("Setting simulation parameters...\n")
sim_params <- list(
  n = 200,      # Sample size
  p = 1500,      # Total number of predictors
  p0 = 20,      # Number of true predictors
  rho = 0.5,    # Correlation for AR-1 covariance matrix
  sigma = 1,    # Standard deviation of the error term
  k = 200,       # Number of simulation replicates
  seed = 2023,  # Seed for reproducibility
  cores = detectCores()    # Number of CPU cores for parallel processing
)

run_config <- list(
  taus = c(0.5, 0.8),
  lambda_grid = seq(0.01, 0.2, length.out = 20),
  zeta_grid = seq(0, 1, length.out = 11)
)


# --- 3. DATA GENERATION ----
cat("Generating simulated data...\n")
set.seed(sim_params$seed)
true_beta <- rep(0, sim_params$p)
true_beta[1:sim_params$p0] <- 1
cov_matrix <- outer(1:sim_params$p, 1:sim_params$p, function(i, j) sim_params$rho^abs(i-j))
X_data <- mvrnorm(n = sim_params$n, mu = rep(0, sim_params$p), Sigma = cov_matrix)
#normal
Y_data <- t(replicate(sim_params$k, as.vector(X_data %*% true_beta +
                                                rnorm(sim_params$n, sd = sim_params$sigma))))
#t3
#Y_data <- t(replicate(sim_params$k, as.vector(X_data %*% true_beta +
#                                                rt(sim_params$n,df = 3))))
cat("Data generation complete.\n")


# --- 4. DEFINE PRIOR SETS ----
cat("Defining prior information sets...\n")
priors <- list(
  S1 = rep(1, sim_params$p), S2 = rep(1, sim_params$p), S3 = rep(1, sim_params$p),
  S4 = rep(1, sim_params$p), S6 = rep(1, sim_params$p)
)
priors$S1[1:20] <- 0                                     # 20 correct variables
priors$S2[c(1:10, (sim_params$p - 1):sim_params$p)] <- 0 # 10 correct, 2 wrong
priors$S6[c(1:10, (sim_params$p - 4):sim_params$p)] <- 0 # 10 correct, 5 wrong
priors$S4[(sim_params$p - 19):sim_params$p] <- 0         # 20 wrong, last 20
priors$S3[16:35] <- 0                                    # 5 correct, 15 wrong


# --- 5. RUN SIMULATIONS ----
cat("Starting simulations...\n")
all_results <- list()

for (tau_val in run_config$taus) {
  cat(sprintf("\n--- Running for tau = %.1f ---\n", tau_val))
  tau_id <- paste0("tau_", tau_val)
  all_results[[tau_id]] <- list()

  # 1. Traditional Penalized QR
  cat("Running Traditional Penalized QR...\n")
  all_results[[tau_id]]$trad_qr <- do.call(
    rbind,
    mclapply(1:sim_params$k, FUN = run_pen_qr_replicate,
             X_input = X_data, Y_input = Y_data, tau = tau_val, mc.cores = sim_params$cores)
  )

  # 2. Loop through each prior set
  for (prior_name in names(priors)) {
    cat(sprintf("Running models for prior set: %s\n", prior_name))

    # 2a. Prior-Informed Penalized QR
    all_results[[tau_id]]$prior_qr[[prior_name]] <- do.call(
      rbind,
      mclapply(1:sim_params$k, FUN = run_pen_qr_replicate,
               X_input = X_data, Y_input = Y_data, tau = tau_val,
               penalty_factor = priors[[prior_name]], mc.cores = sim_params$cores)
    )

    # 2b. Prior-Informed KIQR (using the new API function)
    prior_beta_for_kiqr <- all_results[[tau_id]]$prior_qr[[prior_name]]
    pb_kiqr <- progress_bar$new(total = sim_params$k, format = "Replicate [:bar] :percent")

    kiqr_raw_results <- lapply(1:sim_params$k, function(i) {
      pb_kiqr$tick()
      # Call the new, clean API function for the 2D search.
      KIQR_2D_search(
        X = X_data,
        Y = Y_data[i, ],
        beta_int = rep(0.1, ncol(X_data) + 1),
        beta_p = prior_beta_for_kiqr[i, ],
        tau = tau_val,
        lambda_grid = run_config$lambda_grid,
        zeta_grid = run_config$zeta_grid,
        cores = sim_params$cores
      )
    })

    # Unpack results
    all_results[[tau_id]]$kiqr[[prior_name]] <- list(
      betas = do.call(rbind, lapply(kiqr_raw_results, `[[`, "beta")),
      lambdas = sapply(kiqr_raw_results, `[[`, "lambda"),
      zetas = sapply(kiqr_raw_results, `[[`, "zeta")
    )
  }
}

# --- 6. SAVE RESULTS ----
cat("\nSaving results...\n")
save(
  all_results,
  sim_params,
  true_beta,
  file = "simulation_results_final_n200p1500_norm_k200.Rdata"
)
cat("Results saved.\n")

# --- 1. SETUP ---
cat("Loading required packages and simulation results...\n")
# Make sure all necessary packages are installed
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("knitr")) install.packages("knitr")
if (!require("kableExtra")) install.packages("kableExtra")

library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

# Load the simulation results file
# Note: Ensure the file "simulation_results_final_n200p1500_X1e_k200.Rdata" is in your working directory
load("simulation_results_final_n200p1500_X1e_k200.Rdata")

# --- 2. DEFINE PARAMETERS FOR THE TABLE ---
cat("Defining parameters for the new summary table...\n")

# Scenarios to include in the table
taus_to_report <- c(0.5, 0.8)
priors_to_report <- c("None", "S1", "S2", "S4", "S6")
models_to_report <- c("Trad-QR", "Prior-QR", "KIQR")

# Define true variables, excluding X1 for TP/FP calculation
# Original true variables are 1:20. X1 is index 1.
true_vars_no_x1 <- 2:sim_params$p0

# Define the threshold for considering a coefficient non-zero
selection_threshold <- 1e-4
summary_list <- list()

# --- 3. CALCULATE METRICS FOR EACH SCENARIO ---
cat("Calculating metrics (TP, FP, F1, MSE excluding X1, and X1 Selection Rate)...\n")

for (tau_val in taus_to_report) {

  # Define the true beta vector for this specific tau
  true_beta_for_tau <- true_beta
  true_beta_for_tau[1] <- true_beta[1] * qnorm(tau_val, sd = sim_params$sigma)

  tau_id <- paste0("tau_", tau_val)

  for (prior_name in priors_to_report) {

    # Determine which models to run for this prior
    current_models <- models_to_report
    if (prior_name != "None") {
      current_models <- setdiff(current_models, "Trad-QR")
    } else {
      current_models <- "Trad-QR"
    }

    for (model_name in current_models) {

      # Get the matrix of estimated beta coefficients for the current scenario
      if (model_name == "Trad-QR") {
        beta_matrix <- all_results[[tau_id]]$trad_qr
      } else if (model_name == "Prior-QR") {
        beta_matrix <- all_results[[tau_id]]$prior_qr[[prior_name]]
      } else if (model_name == "KIQR") {
        beta_matrix <- all_results[[tau_id]]$kiqr[[prior_name]]$betas
      }

      # --- Metric Calculations ---

      # 1. X1 Selection Rate (%)
      x1_coeffs <- beta_matrix[, 2] # Column 2 corresponds to X1
      x1_selection_rate <- mean(abs(x1_coeffs) > selection_threshold) * 100

      # 2. TP, FP, and F1 (excluding X1)
      betas_minus_intercept_and_x1 <- beta_matrix[, -(1:2)]

      # NEW: Modified to calculate F1 score in addition to TP and FP
      tp_fp_f1_results <- apply(betas_minus_intercept_and_x1, 1, function(est_coeffs) {
        selected_vars <- which(abs(est_coeffs) > selection_threshold) + 1

        tp <- sum(selected_vars %in% true_vars_no_x1)
        fp <- length(selected_vars) - tp

        # Calculate Precision and Recall for F1 score
        precision <- ifelse((tp + fp) > 0, tp / (tp + fp), 0)
        recall <- tp / length(true_vars_no_x1)

        f1 <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)

        c(TP = tp, FP = fp, F1 = f1) # Return all three metrics
      })

      # NEW: Average F1 score across replicates
      avg_tp <- mean(tp_fp_f1_results["TP", ])
      avg_fp <- mean(tp_fp_f1_results["FP", ])
      avg_f1 <- mean(tp_fp_f1_results["F1", ])

      # 3. MSE (excluding X1)
      true_betas_no_x1 <- true_beta_for_tau[-1]
      sq_errors <- sweep(beta_matrix[, -(1:2)], 2, true_betas_no_x1, "-")^2
      avg_mse <- mean(sq_errors)

      # Append results to the list
      # NEW: Added F1 to the data.frame
      summary_list[[length(summary_list) + 1]] <- data.frame(
        Tau = tau_val,
        Prior = prior_name,
        Model = model_name,
        TP = avg_tp,
        FP = avg_fp,
        F1 = avg_f1,
        MSE = avg_mse,
        X1_Rate = x1_selection_rate
      )
    }
  }
}

# Combine all results into a final data frame
final_summary_df <- bind_rows(summary_list) %>%
  mutate(Prior = factor(Prior, levels = priors_to_report)) %>%
  arrange(Tau, Prior, Model)

# --- 4. GENERATE LATEX TABLE ---
cat("\n\n--- R Code to Generate LaTeX Table ---\n\n")

# Prepare the data frame for kable (formatting and renaming)
# NEW: Formatted and selected the new F1 column
table_df <- final_summary_df %>%
  mutate(
    TP = sprintf("%.2f", TP),
    FP = sprintf("%.2f", FP),
    F1 = sprintf("%.3f", F1), # Format F1 to 3 decimal places
    MSE = sprintf("%.4f", MSE),
    `X1 Rate (%)` = sprintf("%.1f", X1_Rate)
  ) %>%
  select(Tau, Prior, Model, TP, FP, F1, MSE, `X1 Rate (%)`)

# Print the R command to generate the LaTeX code
# NEW: Updated the kable code to include F1
cat('
# This R code generates the LaTeX for the summary table.
# You can copy the output directly into your .tex file.
library(knitr)
library(kableExtra)

# The table data frame is `table_df`
kable(
  table_df,
  format = "latex",
  booktabs = TRUE,
  longtable = FALSE,
  caption = "Simulation Results for Quantile Models (Metrics exclude predictor $X_1$).",
  label = "tab:sim-results",
  col.names = c("Tau", "Prior", "Model", "TP", "FP", "F1", "MSE", "$X_1$ Sel. (\\\\%)"),
  align = "llcccccc",
  escape = FALSE
) %>%
  kable_styling(latex_options = c("striped", "scale_down")) %>%
  collapse_rows(columns = 1:2, valign = "top") %>%
  add_header_above(c(" " = 3, "Performance (excluding $X_1$)" = 4, " " = 1))
')


# --- 5. EXECUTE THE CODE TO SHOW LATEX OUTPUT ---
cat("\n\n--- Generated LaTeX Code ---\n\n")

# Execute the code above to print the LaTeX table
# NEW: Updated the kable code to include F1 in the final output
latex_code <- kable(
  table_df,
  format = "latex",
  booktabs = TRUE,
  longtable = FALSE,
  caption = "Simulation Results for Quantile Models (Metrics exclude predictor $X_1$).",
  label = "tab:sim-results",
  col.names = c("Tau", "Prior", "Model", "TP", "FP", "F1", "MSE", "$X_1$ Sel. (\\%)"),
  align = "llcccccc", # Added a 'c' for the new column
  escape = FALSE
) %>%
  kable_styling(latex_options = c("striped", "scale_down")) %>%
  collapse_rows(columns = 1:2, valign = "top") %>%
  add_header_above(c(" " = 3, "Performance (excluding $X_1$)" = 4, " " = 1)) # Changed span from 3 to 4
table_df
cat(latex_code)

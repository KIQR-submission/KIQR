# SETUP ----
cat("Loading required packages and simulation results...\n")
library(dplyr)
library(tidyr)
library(ggplot2)

# Load the single, consolidated results file
load("simulation_results_final_n200p1500_t3_k200.Rdata")

# -------------------------------------------------------------------------

# HELPER FUNCTION ----
#' Calculates performance metrics for a single simulation run.
compile_metrics <- function(beta_est, true_vars_indices, threshold = 1e-4) {
  # Exclude intercept for variable selection
  selected_vars <- which(abs(beta_est[-1]) > threshold)

  TP <- sum(selected_vars %in% true_vars_indices)
  FP <- length(selected_vars) - TP
  FN <- length(true_vars_indices) - TP

  # F1 Score calculation, handle division by zero
  f1 <- if ((TP + FP + FN) == 0) 0 else (2 * TP) / (2 * TP + FP + FN)

  return(data.frame(TP = TP, FP = FP, F1 = f1))
}

# -------------------------------------------------------------------------

# PROCESS RESULTS ----
cat("Processing results into a tidy data frame...\n")
results_list <- list()
true_variables <- 1:sim_params$p0

# Loop through the nested list of results and build a single data frame
for (tau_id in names(all_results)) {
  tau_val <- as.numeric(gsub("tau_", "", tau_id))

  # Process traditional QR results
  res_trad <- apply(all_results[[tau_id]]$trad_qr, 1, compile_metrics, true_vars_indices = true_variables)
  results_list[[length(results_list) + 1]] <- bind_rows(res_trad, .id = "Replicate") %>%
    mutate(Tau = tau_val, Model = "Trad-QR", Prior = "None")

  # Process prior QR and KIQR results
  for (prior_name in names(all_results[[tau_id]]$prior_qr)) {
    # Prior QR
    res_prior <- apply(all_results[[tau_id]]$prior_qr[[prior_name]], 1, compile_metrics, true_vars_indices = true_variables)
    results_list[[length(results_list) + 1]] <- bind_rows(res_prior, .id = "Replicate") %>%
      mutate(Tau = tau_val, Model = "Prior-QR", Prior = prior_name)

    # KIQR - Access the 'betas' matrix inside the list
    res_kiqr <- apply(all_results[[tau_id]]$kiqr[[prior_name]]$betas, 1, compile_metrics, true_vars_indices = true_variables)
    results_list[[length(results_list) + 1]] <- bind_rows(res_kiqr, .id = "Replicate") %>%
      mutate(Tau = tau_val, Model = "KIQR", Prior = prior_name)
  }
}

# Combine all data frames into one
final_results_df <- bind_rows(results_list) %>%
  mutate(Replicate = as.integer(Replicate),
         Prior = factor(Prior, levels = c("None", "S1", "S2", "S6", "S3", "S4", "S5" )))

# -------------------------------------------------------------------------

# VISUALIZE RESULTS ----
cat("Generating and saving plots...\n")

# Filter data to match original plot (Priors S1-S4)
plot_df <- final_results_df %>%
  filter(Prior %in% c("None", "S1", "S2", "S6", "S4")) %>%
  filter(Tau %in% c(0.5, 0.8)) %>%
  mutate(
    Prior = case_when(
      Prior == "S6" ~ "S3",
      TRUE ~ Prior # This line keeps all other values as they are
    )
  )

library(forcats)

# 1. Reorder the factor levels in the plot_df
plot_df$Model <- fct_relevel(plot_df$Model,
                             "Trad-QR", "Prior-QR", "KIQR")


# Create a combined plot with facets for each tau
f1_boxplot_grayscale <- ggplot(plot_df, aes(x = Prior, y = F1, fill = Model)) + # Keep fill for box color
  geom_boxplot() +
  # Add geom_point for mean/median with shape for better differentiation
  stat_summary(fun = median, geom = "point", aes(shape = Model), # Use median for consistency with boxplot
               size = 3, position = position_dodge(width = 0.75)) + # Dodge points to match boxplot positions
  facet_wrap(~ paste("Tau =", Tau), ncol = 2) +
  # Use a grayscale palette for fills
  scale_fill_grey(start = 0.8, end = 0.2) + # Adjust start/end for light to dark gray
  # Manually define shapes for better control (up to 6 shapes are easily distinguishable)
  scale_shape_manual(values = c("Trad-QR" = 24, "Prior-QR" = 22, "KIQR" = 21)) + # Filled shapes
  scale_fill_manual(values = c("Trad-QR" = "#F8766D", "Prior-QR" = "#00BA38", "KIQR" = "#619CFF")) +
  labs(
    x = "Prior Information Set",
    y = "F1 Score",
    fill = "Model Type",
    shape = "Model Type" # Add shape legend title
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey80"), # Light grey background for facets
    legend.position = "bottom",
    # Further customization for a "journal-ready" look
    panel.grid.major = element_line(linetype = "dotted", color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.7)
  )

# Print the plot to the screen
print(f1_boxplot_grayscale)

# Save the plot to a file
ggsave("F1_n200p1500_t3_bwc.png", plot = f1_boxplot_grayscale, width = 12, height = 7, dpi = 300)

cat("Analysis complete. Plot saved to 'F1_Scores_Comparison.png'.\n")

# Generate a summary table
summary_table <- final_results_df %>% group_by(Tau, Prior, Model) %>% summarise(TP = mean(TP),
                                                               FP = mean(FP),
                                                               F1 = mean(F1))


# If hete X, calculate X1 pick up rate:
cat("\nCalculating selection percentage for X1...\n")

# Define the threshold for selection, consistent with compile_metrics
selection_threshold <- 1e-4
x1_selection_list <- list()

# Loop through the results structure to check for X1 selection
for (tau_id in names(all_results)) {
  tau_val <- as.numeric(gsub("tau_", "", tau_id))

  # 1. Traditional QR
  # The coefficient for X1 is in the second column (first is the intercept)
  betas_x1_trad <- all_results[[tau_id]]$trad_qr[, 2]
  selection_pct_trad <- mean(abs(betas_x1_trad) > selection_threshold) * 100
  x1_selection_list[[length(x1_selection_list) + 1]] <- data.frame(
    Tau = tau_val, Model = "Trad-QR", Prior = "None", SelectionPct = selection_pct_trad
  )

  # 2. Prior QR and KIQR
  for (prior_name in names(all_results[[tau_id]]$prior_qr)) {
    # Prior-QR
    betas_x1_prior <- all_results[[tau_id]]$prior_qr[[prior_name]][, 2]
    selection_pct_prior <- mean(abs(betas_x1_prior) > selection_threshold) * 100
    x1_selection_list[[length(x1_selection_list) + 1]] <- data.frame(
      Tau = tau_val, Model = "Prior-QR", Prior = prior_name, SelectionPct = selection_pct_prior
    )

    # KIQR
    betas_x1_kiqr <- all_results[[tau_id]]$kiqr[[prior_name]]$betas[, 2]
    selection_pct_kiqr <- mean(abs(betas_x1_kiqr) > selection_threshold) * 100
    x1_selection_list[[length(x1_selection_list) + 1]] <- data.frame(
      Tau = tau_val, Model = "KIQR", Prior = prior_name, SelectionPct = selection_pct_kiqr
    )
  }
}

# Combine into a single data frame
x1_selection_df <- bind_rows(x1_selection_list)

# Reshape the data for a clean table output, rounding the percentages
x1_summary_table <- x1_selection_df %>%
  mutate(SelectionPct = round(SelectionPct, 1)) %>% # Round for cleaner output
  pivot_wider(
    names_from = Model,
    values_from = SelectionPct
  ) %>%
  arrange(Tau, Prior) %>%
  select(Tau, Prior, `Trad-QR`, `Prior-QR`, KIQR) # Reorder columns

# Print the final table using kable for nice formatting
cat("\n--- Percentage of Replicates where X1 was Selected (%) ---\n")
print(knitr::kable(x1_summary_table, format = "pipe", caption = "Selection Percentage for X1"))

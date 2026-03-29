#' Fit a KIQR Model using a Proximal Gradient Algorithm
#'
#' This is the primary function for fitting the Knowledge-Informed Quantile Regression
#' (KIQR) model using a Proximal Gradient Descent (PGD) algorithm. It is generally
#' preferred for its stability and handling of the non-differentiable SCAD penalty.
#'
#' @param X The design matrix (predictors), without an intercept column.
#' @param Y The primary response vector.
#' @param Y_p The auxiliary (prior knowledge) response vector.
#' @param beta_int Optional initial values for the beta coefficients. If NULL, they are generated automatically.
#' @param beta_p Prior beta coefficients from the auxiliary model.
#' @param zeta The transfer learning parameter (0 <= zeta <= 1). zeta=0 corresponds to standard QR,
#'   zeta=1 corresponds to using only prior information.
#' @param gamma The smoothing parameter for the Huber loss function.
#' @param lambda The regularization parameter for the SCAD penalty.
#' @param tau The quantile to be estimated.
#' @param a The second tuning parameter for the SCAD penalty (default is 3.7).
#' @param max.iter The maximum number of iterations.
#' @param tol The convergence tolerance.
#' @param penalty_factor A vector of penalty multipliers for each coefficient.
#' @param prior_pen_factor Penalty factors for the prior-only model (when zeta > 1).
#' @param half_life The half-life for the decaying step size.
#' @return A numeric vector of the estimated beta coefficients.
#' @export
KIQR_R_iteration_zeta_proximal <- function(X, Y, Y_p, beta_int = NULL, beta_p = NULL,
                                           zeta, gamma, lambda, tau, a = 3.7,
                                           max_iter = 1000, tol = 1e-4, penalty_factor = NULL,
                                           prior_pen_factor = NULL, half_life = 100) {
  p <- NCOL(X)
  X_design <- cbind(1, X)

  # --- Initialize Beta Coefficients ---
  # check beta_int
  if(is.null(beta_int)){
    fit_y <- rq.pen(x = X,y = as.vector(Y),tau = tau,penalty = "SCAD")
    beta_int_y <- as.vector(coef(fit_y, lambda = select_lambda_hbic(fit_y)))
  }else{
    beta_int_y <- beta_int
  }
  # check beta_p
  if(is.null(beta_p)){
    fit_p <- rq.pen(x = X,y = as.vector(Y_p),tau = tau)
    beta_int_yp <- as.vector(coef(fit_p, lambda = select_lambda_hbic(fit_p)))
  }else{
    beta_p <- beta_p
  }
  # weighted beta_int
  # browser()
  beta_int = zeta*beta_p + (1-zeta)*beta_int_y

  # --- Setup Penalty Factors ---
  if (is.null(penalty_factor)) {
    penalty_factor <- c(0, rep(1, p))
  } else {
    penalty_factor <- c(0, penalty_factor)
  }
  if (is.null(prior_pen_factor)) {
    prior_pen_factor <- c(0, rep(1, p))
  } else {
    prior_pen_factor <- c(0, prior_pen_factor)
  }

  # --- Call Solver ---
  if (zeta > 1) { # a test case, to be deleted
    beta_res <- solver_pgd(beta = beta_int, Y = Y, X = X_design, Y_p = Y_p,
                           tau = tau, zeta = 0, lambda = lambda,
                           max_iter = max_iter, gamma = gamma, tol = tol, a = a,
                           penalty_factor = prior_pen_factor, halflife = half_life)
  } else {
    beta_res <- solver_pgd(beta = beta_int, Y = Y, X = X_design, Y_p = Y_p,
                           tau = tau, zeta = zeta, lambda = lambda,
                           max_iter = max_iter, gamma = gamma, tol = tol, a = a,
                           penalty_factor = penalty_factor, halflife = half_life)
  }
  return(beta_res)
}


#' Fit a KIQR Model with 2D Grid Search Tuning
#'
#' Performs a comprehensive 2D grid search over lambda and zeta to find the
#' optimal KIQR model, evaluated by HBIC. This is the recommended high-level
#' function for most use cases.
#'
#' @param X The predictor matrix.
#' @param Y The response vector for the current replicate.
#' @param beta_p The prior beta coefficient vector for this replicate.
#' @param tau The quantile level.
#' @param lambda_grid A numeric vector of lambda values to test.
#' @param zeta_grid A numeric vector of zeta values to test.
#' @param penalty_factor Optional penalty factor vector.
#' @param cores The number of CPU cores to use for parallelizing the grid search.
#' @param gamma The smoothing parameter for the Huber loss function.
#' @param half_life The half-life for the decaying step size.
#' @return A list containing the best beta vector (`beta`), the optimal lambda (`lambda`),
#'   and the optimal zeta (`zeta`).
#' @export
KIQR_2D_search <- function(X, Y, beta_p, tau, lambda_grid, zeta_grid,
                           beta_int = NULL,
                           penalty_factor = NULL, cores = 1,
                           gamma = 0.01, half_life = 100, max_iter = 1000) {

  # 1. Create the 2D parameter grid
  param_grid <- expand.grid(lambda = lambda_grid, zeta = zeta_grid)

  # Generate the prior-based response
  Y_p <- cbind(1, X) %*% beta_p

  # 2. Parallelize the grid search
  beta_candidates_list <- parallel::mclapply(1:nrow(param_grid), function(j) {
    params <- param_grid[j, ]
    KIQR_R_iteration_zeta_proximal(
      beta_int = beta_int,
      X = X,
      Y = Y,
      Y_p = Y_p,
      gamma = gamma,
      lambda = params$lambda,
      zeta = params$zeta,
      tau = tau,
      half_life = half_life,
      beta_p = beta_p,
      penalty_factor = penalty_factor,
      max_iter = max_iter
    )
  }, mc.cores = cores)

  beta_candidates <- do.call(rbind, beta_candidates_list)

  # 3. Find the best combination using HBIC
  hbic_scores <- hbic_from_beta(X, Y, beta_candidates, tau)
  best_index <- which.min(hbic_scores)

  best_beta <- beta_candidates[best_index, ]
  best_params <- param_grid[best_index, ]

  # 4. Return results as a list
  return(list(
    beta = best_beta,
    lambda = best_params$lambda,
    zeta = best_params$zeta
  ))
}

#' Fit a KIQR Model with 2D k-Fold Cross-Validation Tuning (Parallel over Params)
#'
#' Performs a comprehensive 2D grid search over lambda and zeta to find the
#' optimal KIQR model, evaluated by k-fold cross-validation. This version
#' parallelizes the search over the (lambda, zeta) parameter grid, which is
#' efficient when the number of parameter combinations is large.
#'
#' @param X The predictor matrix.
#' @param Y The response vector.
#' @param beta_p The prior beta coefficient vector.
#' @param tau The quantile level.
#' @param lambda_grid A numeric vector of lambda values to test.
#' @param zeta_grid A numeric vector of zeta values to test.
#' @param nfolds The number of folds for cross-validation. Default is 5.
#' @param beta_int Optional initial beta vector for the optimization.
#' @param penalty_factor Optional penalty factor vector.
#' @param cores The number of CPU cores to use for parallelizing the grid search
#'   (parallelizes over the parameter combinations).
#' @param gamma The smoothing parameter for the Huber loss function.
#' @param half_life The half-life for the decaying step size.
#' @param max_iter The maximum number of iterations for the proximal gradient descent.
#' @return A list containing:
#'   \item{beta.min}{The beta vector from the model fit with `lambda.min` and `zeta.min`.}
#'   \item{lambda.min}{The optimal lambda (minimizing mean CV loss).}
#'   \item{zeta.min}{The optimal zeta (minimizing mean CV loss).}
#'   \item{beta.1se}{The beta vector from the model fit with `lambda.1se` and `zeta.1se`.}
#'   \item{lambda.1se}{The lambda selected by the 1-SE rule (largest lambda within 1 SE of the minimum).}
#'   \item{zeta.1se}{The zeta corresponding to `lambda.1se`.}
#'   \item{param_grid}{A data frame of the grid search results, including `cv_mean` and `cv_se` for each parameter pair.}
#'   \item{nfolds}{The number of folds used.}
#' @export
KIQR_2D_CV <- function(X, Y, beta_p, tau,lambda_grid, zeta_grid,
                       nfolds = 5,
                       beta_int = NULL,
                       penalty_factor = NULL, cores = 1,
                       gamma = 0.01, half_life = 100, max_iter = 1000) {

  # --- 1. Helper Function for Quantile Check Loss ---
  quantile_check_loss <- function(y, y_hat, tau) {
    residual <- y - y_hat
    loss <- ifelse(residual >= 0, tau * residual, (tau - 1) * residual)
    return(mean(loss))
  }

  # --- 2. Setup Folds and Parameter Grid ---
  n <- length(Y)
  if (n != nrow(X)) {
    stop("Length of Y must match number of rows in X.")
  }

  # Create fold indices
  fold_indices <- sample(cut(seq(1, n), breaks = nfolds, labels = FALSE))

  # Create the 2D parameter grid
  param_grid <- expand.grid(lambda = lambda_grid, zeta = zeta_grid)
  n_params <- nrow(param_grid)

  # --- 3. Run k-Fold Cross-Validation (Parallel over parameter grid) ---
  # We parallelize over the (lambda, zeta) combinations.
  # Each core takes a combo and runs the full k-fold CV for it.
  cv_results_list <- parallel::mclapply(1:n_params, function(j) {

    params <- param_grid[j, ]

    # Vector to store losses for this parameter combo (one per fold)
    fold_losses <- numeric(nfolds)

    # Loop over every fold
    for (i in 1:nfolds) {
      # Define training and validation sets for this fold
      val_idx <- which(fold_indices == i)
      train_idx <- which(fold_indices != i)

      X_train <- X[train_idx, , drop = FALSE]
      Y_train <- Y[train_idx]
      X_val <- X[val_idx, , drop = FALSE]
      Y_val <- Y[val_idx]

      # 1. Fit model on training data
      # Calculate the prior-based response *for the training set*
      Y_p_train <- (cbind(1, X_train) %*% beta_p)

      beta_hat <- tryCatch({
        KIQR_R_iteration_zeta_proximal(
          beta_int = beta_int, # Pass the same initial guess
          X = X_train,
          Y = Y_train,
          Y_p = Y_p_train,
          gamma = gamma,
          lambda = params$lambda,
          zeta = params$zeta,
          tau = tau,
          half_life = half_life,
          beta_p = beta_p, # beta_p is (p+1)-dim, pass as is
          penalty_factor = penalty_factor,
          max_iter = max_iter
        )
      }, error = function(e) {
        # In case of convergence failure, return NA
        return(NA)
      })

      if (any(is.na(beta_hat))) {
        fold_losses[i] <- NA
        next
      }

      # 2. Predict on validation data
      Y_hat_val <- cbind(1, X_val) %*% beta_hat

      # 3. Calculate and store quantile loss
      fold_losses[i] <- quantile_check_loss(Y_val, Y_hat_val, tau)
    }

    # After all folds are done, calculate mean and sd for this param combo
    n_valid_folds <- sum(!is.na(fold_losses))
    mean_loss <- mean(fold_losses, na.rm = TRUE)
    sd_loss <- sd(fold_losses, na.rm = TRUE)

    # Return the aggregated results for this (lambda, zeta) pair
    return(list(
      fold_losses = fold_losses,
      mean_loss = mean_loss,
      sd_loss = sd_loss,
      n_valid = n_valid_folds
    ))

  }, mc.cores = cores)

  # --- 4. Process CV Results ---
  # Extract the results for all parameter combinations
  fold_losses <- sapply(cv_results_list, function(x) x$fold_losses)
  cv_mean_loss <- sapply(cv_results_list, function(x) x$mean_loss)
  cv_sd_loss <- sapply(cv_results_list, function(x) x$sd_loss)
  n_valid <- sapply(cv_results_list, function(x) x$n_valid)

  # Calculate standard error of the mean
  cv_std_err <- cv_sd_loss / sqrt(n_valid)

  # Handle cases with 0 or 1 valid folds (where sd is NA)
  cv_std_err[is.na(cv_std_err)] <- NA

  # Store results in the parameter grid
  param_grid$cv_mean <- cv_mean_loss
  param_grid$cv_se <- cv_std_err

  # --- 5. Find Best Parameters (.min) ---
  best_index_min <- which.min(param_grid$cv_mean)
  if (length(best_index_min) == 0 || all(is.na(param_grid$cv_mean))) {
    stop("All CV fits failed or produced NAs. Check parameters or model stability.")
  }

  best_params_min <- param_grid[best_index_min, ]
  lambda_min <- best_params_min$lambda
  zeta_min <- best_params_min$zeta
  min_cv_loss <- best_params_min$cv_mean

  # --- 6. Find Best Parameters (.1se) ---
  # 1-SE Rule: Find the most "parsimonious" model whose mean CV loss is
  # within one standard error of the minimum.
  # "Parsimonious" = largest lambda, then largest zeta.

  # Ensure threshold is valid even if se is NA (e.g., 1 fold)
  se_min <- param_grid$cv_se[best_index_min]
  if (is.na(se_min)) { se_min <- 0 }

  threshold <- min_cv_loss + se_min

  # Get all candidates within the 1-SE threshold
  param_grid_1se_candidates <- param_grid[param_grid$cv_mean <= threshold & !is.na(param_grid$cv_mean), ]

  if (nrow(param_grid_1se_candidates) > 0) {
    # Find the largest lambda among them
    max_lambda_1se <- max(param_grid_1se_candidates$lambda)

    # Get all candidates with this largest lambda
    final_1se_candidates <- subset(param_grid_1se_candidates, lambda == max_lambda_1se)

    # From these, pick the one with the largest zeta
    best_params_1se <- final_1se_candidates[which.max(final_1se_candidates$zeta), ]

    lambda_1se <- best_params_1se$lambda
    zeta_1se <- best_params_1se$zeta
  } else {
    # Fallback in case 1-SE rule fails
    lambda_1se <- lambda_min
    zeta_1se <- zeta_min
  }

  # --- 7. Refit Models on Full Data ---

  # Generate the full prior-based response
  Y_p_full <- cbind(1, X) %*% beta_p

  # Fit the final .min model
  beta_min <- KIQR_R_iteration_zeta_proximal(
    beta_int = beta_int,
    X = X,
    Y = Y,
    Y_p = Y_p_full,
    gamma = gamma,
    lambda = lambda_min,
    zeta = zeta_min,
    tau = tau,
    half_life = half_life,
    beta_p = beta_p,
    penalty_factor = penalty_factor,
    max_iter = max_iter
  )

  # Fit the final .1se model
  beta_1se <- KIQR_R_iteration_zeta_proximal(
    beta_int = beta_int,
    X = X,
    Y = Y,
    Y_p = Y_p_full,
    gamma = gamma,
    lambda = lambda_1se,
    zeta = zeta_1se,
    tau = tau,
    half_life = half_life,
    beta_p = beta_p,
    penalty_factor = penalty_factor,
    max_iter = max_iter
  )

  # --- 8. Return Results ---
  return(list(
    beta.min = beta_min,
    lambda.min = lambda_min,
    zeta.min = zeta_min,

    beta.1se = beta_1se,
    lambda.1se = lambda_1se,
    zeta.1se = zeta_1se,

    param_grid = param_grid,
    fold_losses = fold_losses,
    nfolds = nfolds
  ))
}

#' Fit a KIQR Model using a Gradient Descent Algorithm
#'
#' Fits the KIQR model using a standard gradient descent approach with hard
#' thresholding for variable selection. The PGD version (`KIQR_R_iteration_zeta_proximal`)
#' is often preferred.
#'
#' @inheritParams KIQR_R_iteration_zeta_proximal
#' @param cut_off A hard threshold below which coefficients are set to zero after each iteration.
#' @return A numeric vector of the estimated beta coefficients.
#' @export
KIQR_R_iteration_zeta <- function(X, Y, Y_p, beta_int = NULL, beta_p = NULL,
                                  zeta, gamma, lambda, tau, a = 3.7,
                                  max.iter = 10000, tol = 1e-4,
                                  cut_off = 1e-5, penalty_factor = NULL,
                                  prior_pen_factor = NULL, half_life = 100) {
  p <- NCOL(X)
  X_design <- cbind(1, X)

  # --- Initialize Beta Coefficients ---
  if (is.null(beta_int)) {
    fit_0 <- rqPen::rq.pen(x = X, y = as.vector(Y), tau = tau)
    beta_int_y <- as.vector(stats::coef(fit_0, lambda = select_lambda_hbic(fit_0)))
    if (!is.null(beta_p)) {
      beta_int <- zeta * beta_p + (1 - zeta) * beta_int_y
    } else {
      fit_1 <- rqPen::rq.pen(x = X, y = as.vector(Y_p), tau = tau)
      beta_int_yp <- as.vector(stats::coef(fit_1, lambda = select_lambda_hbic(fit_1)))
      beta_int <- zeta * beta_int_yp + (1 - zeta) * beta_int_y
    }
  }

  # --- Setup Penalty Factors ---
  if (is.null(penalty_factor)) penalty_factor <- c(0, rep(1, p))
  else penalty_factor <- c(0, penalty_factor)

  if (is.null(prior_pen_factor)) prior_pen_factor <- c(0, rep(1, p))
  else prior_pen_factor <- c(0, prior_pen_factor)


  # --- Call Solver ---
  if (zeta > 1) {
    beta_res <- solver_gd(beta = beta_int, Y = Y, X = X_design, Y_p = Y_p,
                          tau = tau, zeta = 0, lambda = lambda, max_iter = max.iter,
                          gamma = gamma, tol = tol, cut_off = cut_off,
                          penalty_factor = prior_pen_factor, halflife = half_life)
  } else {
    beta_res <- solver_gd(beta = beta_int, Y = Y, X = X_design, Y_p = Y_p,
                          tau = tau, zeta = zeta, lambda = lambda, max_iter = max.iter,
                          gamma = gamma, tol = tol, cut_off = cut_off,
                          penalty_factor = penalty_factor, halflife = half_life)
  }
  return(beta_res)
}

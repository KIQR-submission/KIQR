#' Core HBIC Calculation
#'
#' A helper function that computes the Hannan-Quinn-type Information Criterion (HBIC).
#'
#' @param Rho The mean check loss.
#' @param NR The number of observations.
#' @param nvar The number of non-zero variables (excluding the intercept).
#' @param p The total number of predictors.
#' @param n The sample size.
#' @return The calculated HBIC value.
hbic_core <- function(Rho, NR, nvar, p, n) {
  log(Rho * NR) + nvar * log(log(p)) * log(n) / n
}

#' Calculate HBIC from Beta Coefficients
#'
#' Computes the HBIC for a given set of coefficients. This function can handle
#' a single coefficient vector or a matrix where each row is a coefficient vector.
#'
#' @param X The design matrix (without intercept).
#' @param Y The response vector.
#' @param beta A numeric vector or matrix of coefficients (including intercept).
#' @param tau The quantile level.
#' @param threshold The numerical tolerance to determine if a coefficient is non-zero.
#' @return A single HBIC value if `beta` is a vector, or a vector of HBIC values
#'   if `beta` is a matrix.
#' @export
hbic_from_beta <- function(X, Y, beta, tau, threshold = 0) {
  # If beta is a matrix, apply this function to each row
  if (is.matrix(beta)) {
    return(apply(beta, 1, function(beta_row) {
      hbic_from_beta(X = X, Y = Y, beta = beta_row, tau = tau, threshold = threshold)
    }))
  }

  # Logic for a single beta vector
  n <- NROW(X)
  p <- NCOL(X)
  X_design <- cbind(1, X)

  residual <- Y - X_design %*% beta
  rho <- mean(check_loss(u = residual, tau = tau))
  # nvar should not include the intercept. The original code included it.
  nvar <- sum(abs(beta[-1]) > threshold) # change it!
  nvar <- sum(abs(beta) > threshold) # change it!

  hbic_core(Rho = rho, NR = n, nvar = nvar, p = p, n = n)
}

#' Select Lambda via HBIC from an rqPen Fit
#'
#' Selects the optimal lambda from a `rq.pen` model fit object based on the
#' minimum HBIC value.
#'
#' @param qr_fit A model object returned by `rqPen::rq.pen()`.
#' @param threshold An optional threshold to determine non-zero coefficients. If NULL,
#'   the model's own count (`nzero`) is used.
#' @param HBIC.min If TRUE, returns both the selected lambda and the minimum HBIC value.
#' @return The selected lambda, or a vector c(lambda, min_hbic) if `HBIC.min` is TRUE.
#' @export
select_lambda_hbic <- function(qr_fit, threshold = NULL, HBIC.min = NULL) {
  n <- qr_fit$n
  p <- qr_fit$p
  fit_rho <- qr_fit$models[[1]]$rho

  if (is.null(threshold)) {
    n_nonzero <- qr_fit$models[[1]]$nzero
  } else {
    # Exclude intercept column (first column) from count
    coeffs <- qr_fit$models[[1]]$coefficients[-1, , drop = FALSE]
    n_nonzero <- apply(coeffs, 2, function(x) sum(abs(x) > threshold))
  }
  HBIC <- log(fit_rho * n) + n_nonzero * log(log(p)) * log(n) / n
  #print(cbind(qr_fit$lambda,log(fit_rho * n),n_nonzero * log(log(p)) * log(n) / n,HBIC))
  lambda_selected <- qr_fit$lambda[which.min(HBIC)]

  if (!is.null(HBIC.min)) {
    res <- c(lambda_selected = lambda_selected, min_hbic = min(HBIC))
  } else {
    res <- lambda_selected
  }
  return(res)
}

#' Select Lambda via QBIC from an rqPen Fit
#'
#' Selects the optimal lambda based on the minimum QBIC.
#'
#' @param qr_fit A model object returned by `rqPen::rq.pen()`.
#' @return The selected lambda value.
#' @export
select_lambda_qbic <- function(qr_fit) {
  n <- qr_fit$n
  p <- qr_fit$p
  rho_vals <- qr_fit$models[[1]]$rho
  n_zeros <- qr_fit$models[[1]]$nzero

  QBIC <- log(rho_vals * n) + n_zeros * log(p) * log(log(n)) / (2 * n)
  lambda_selected <- qr_fit$lambda[which.min(QBIC)]
  return(lambda_selected)
}

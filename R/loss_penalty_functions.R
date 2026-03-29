#' Calculate Quantile Check Loss
#'
#' Computes the check loss function, which is fundamental for quantile regression.
#' The loss is defined as $\rho_{\tau}(u) = u(\tau - I(u < 0))$.
#'
#' @param u A numeric vector of residuals ($Y - \hat{Y}$).
#' @param tau The quantile level, a single numeric value between 0 and 1.
#' @return A numeric vector of the same length as `u` containing the check loss values.
#' @export
check_loss <- function(u, tau) {
  u * (tau - (u < 0))
}

#' Calculate the Derivative of the Huberized Check Loss
#'
#' A smoothed, vectorized version of the check loss derivative using the Huber method.
#' This is used to provide a smooth gradient for optimization algorithms.
#'
#' @param u A numeric vector of residuals.
#' @param tau The quantile level (0 < tau < 1).
#' @param gamma The Huber smoothing parameter. Values of `u` with absolute value less than `gamma` are treated quadratically.
#' @return A numeric vector containing the gradient values.
#' @export
huber_check_loss_deriv <- function(u, tau, gamma) {
  is_le_gamma <- abs(u) <= gamma
  grad <- (2 * tau - 1 + (u / gamma) * is_le_gamma + sign(u) * (!is_le_gamma)) / 2
  grad[u == 0] <- 0 # Explicitly handle the non-differentiable point at 0
  return(grad)
}

#' Calculate the Derivative of the SCAD Penalty
#'
#' Computes the derivative of the Smoothly Clipped Absolute Deviation (SCAD) penalty.
#' This is a vectorized function.
#'
#' @param u A numeric vector of coefficients.
#' @param lambda The tuning parameter for the penalty.
#' @param a The second tuning parameter of the SCAD penalty (typically 3.7).
#' @return A numeric vector of the same length as `u` containing the SCAD derivative values.
#' @export
qscad_deriv <- function(u, lambda, a = 3.7) {
  abs_u <- abs(u)
  # Using nested ifelse for vectorization
  grad <- ifelse(abs_u <= lambda,
                 lambda,
                 ifelse(abs_u <= a * lambda,
                        pmax(0, a * lambda - abs_u) / (a - 1),
                        0))
  return(grad)
}

#' Proximal Operator for the SCAD Penalty
#'
#' Applies the proximal operator for the SCAD penalty, a key step in proximal
#' gradient descent algorithms.
#'
#' @param z The input vector (typically the result of a gradient step).
#' @param lambda The main penalty tuning parameter.
#' @param zeta The transfer learning tuning parameter.
#' @param a The second SCAD tuning parameter (must be > 2).
#' @param gamma_vec A vector of step sizes corresponding to each element of `z`.
#' @return The result of the proximal operator applied to `z`.
#' @export
scad_prox <- function(z, lambda, zeta, a, gamma_vec) {
  if (a <= 2) {
    stop("SCAD parameter 'a' must be > 2 for the standard proximal formula.")
  }

  u_star <- numeric(length(z))
  abs_z <- abs(z)
  lambda_g <- (1 - zeta) * lambda * gamma_vec
  gamma_vec <- gamma_vec * (1 - zeta)

  # Vectorized conditions
  cond1 <- abs_z <= lambda + lambda_g
  cond2 <- !cond1 & (abs_z <= a * lambda)
  cond3 <- !cond1 & !cond2

  # Condition 1: Soft-thresholding
  u_star[cond1] <- sign(z[cond1]) * pmax(0, abs_z[cond1] - lambda_g[cond1])

  # Condition 2: More complex threshold
  denominator_term <- 1 - gamma_vec[cond2] / (a - 1)
  if (any(denominator_term <= 1e-8)) {
    warning("SCAD prox condition gamma < a-1 may be violated. Result might be unstable.")
  }
  thresh_val <- gamma_vec[cond2] * a * lambda / (a - 1)
  soft_thresh_interim <- sign(z[cond2]) * pmax(0, abs_z[cond2] - thresh_val)
  u_star[cond2] <- soft_thresh_interim / denominator_term

  # Condition 3: Pass-through
  u_star[cond3] <- z[cond3]

  # Handle cases with no penalty
  u_star[gamma_vec == 0] <- z[gamma_vec == 0]

  return(u_star)
}

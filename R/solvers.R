#' Proximal Gradient Descent (PGD) Solver for KIQR
#'
#' Core iterative algorithm using Proximal Gradient Descent to solve the KIQR objective.
#' Sparsity is induced by the `scad_prox` operator.
#'
#' @inheritParams KIQR_R_iteration_zeta_proximal
#' @param beta Initial beta coefficient vector.
#' @param X The design matrix (with intercept).
#' @param Y The primary response vector.
#' @param Y_p The auxiliary response vector.
#' @return The optimized coefficient vector.
#' @keywords internal
solver_pgd <- function(beta, Y, X, Y_p, tau, gamma, zeta, tol, max_iter, lambda, a,
                       penalty_factor, halflife) {
  n <- nrow(X)
  p <- ncol(X)
  beta_fit <- beta

  # Step size calculation
  xi <- 2 * colSums(X^2) / (n * gamma) #used to be 2 and 30
  xi[xi <= 1e-8] <- 1e-6 # Ensure step sizes are positive and non-zero

  for (i in 1:max_iter) {
    beta_old <- beta_fit
    #browser()
    #print(length(which(beta_old!=0)))
    #print(beta_old[which(beta_old!=0)])
    # Calculate residuals and gradients of the smooth loss part
    Y_fit <- X %*% beta_fit
    residual_obs <- Y - Y_fit
    residual_p <- Y_p - Y_fit
    deriv_obs <- huber_check_loss_deriv(residual_obs, tau, gamma)
    deriv_p <- huber_check_loss_deriv(residual_p, tau, gamma)
    L_obs <- colMeans(as.vector(deriv_obs) * X)
    L_p <- colMeans(as.vector(deriv_p) * X)
    neg_gradient_smooth <- (1 - zeta) * L_obs + zeta * L_p

    # Proximal Gradient Descent Step
    step_size_vec <- (0.5 ^ (i / halflife)) / xi
    z <- beta_fit + step_size_vec * neg_gradient_smooth
    gamma_prox <- step_size_vec * penalty_factor
    beta_fit <- scad_prox(z = z, lambda = lambda, a = a, gamma_vec = gamma_prox, zeta = zeta)

    # Check for convergence
    max_update <- max(abs(beta_fit - beta_old))
    if (any(is.nan(beta_fit)) || any(is.infinite(beta_fit))) {
      warning("NaN or Inf detected in coefficients at iteration ", i, ". Algorithm diverged.")
      return(beta_old) # Return last valid estimate
    }
    if (max_update < tol && i > 4) {
      return(beta_fit)
    }
  }
  message("PGD solver did not converge after ", max_iter, " iterations.")
  return(beta_fit)
}


#' Gradient Descent (GD) Solver for KIQR
#'
#' Core iterative algorithm using a gradient descent variant.
#'
#' @inheritParams KIQR_R_iteration_zeta
#' @param beta Initial beta coefficient vector.
#' @param X The design matrix (with intercept).
#' @param Y The primary response vector.
#' @param Y_p The auxiliary response vector.
#' @return The optimized coefficient vector.
#' @keywords internal
solver_gd <- function(beta, Y, X, Y_p, tau, gamma, zeta, tol, max_iter, lambda,
                      cut_off, penalty_factor, halflife) {
  n <- nrow(X)
  p <- ncol(X)
  beta_fit <- beta

  xi <- 2 * colSums(X^2) / (p * n * gamma)
  xi[xi <= 1e-8] <- 1e-6 # Ensure positivity

  for (i in 1:max_iter) {
    Y_fit <- X %*% beta_fit
    residual_obs <- Y - Y_fit
    residual_p <- Y_p - Y_fit

    deriv_obs <- huber_check_loss_deriv(residual_obs, tau, gamma)
    deriv_p <- huber_check_loss_deriv(residual_p, tau, gamma)
    P_deriv <- qscad_deriv(beta_fit, lambda)

    L_obs <- colMeans(as.vector(deriv_obs) * X)
    L_p <- colMeans(as.vector(deriv_p) * X)

    update <- (0.5 ^ (i/halflife) ) * (1 / xi) * ((1 - zeta) * L_obs +
                                                    zeta * L_p - sign(beta_fit) * penalty_factor * P_deriv)

    beta_fit <- beta_fit + update
    beta_fit <- vec_cut_off(beta_fit, cut_off) # Hard thresholding

    if (max(abs(update)) < tol && i > 4) {
      return(beta_fit)
    }
  }
  message("GD solver did not converge after ", max_iter, " iterations.")
  if(max_iter <=4){
    message("Might because max_iter is set smaller than 4, which is the minimum.")
  }
  return(beta_fit)
}

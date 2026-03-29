###################################################
### Refit the model with the selected variables ###
###################################################

### Inputs:
# x: design matrix
# y: response vector
# beta: original coeffecients
### Output:
# updated coefficients

update.beta <- function(x, y, beta) {
  pos <- which(abs(beta[-1]) > 1e-10)
  beta.refit <- rep(0, length(beta))
  if (length(pos) >= 1) {
    refit <- lm(y ~ x[, pos])
    beta.refit[c(1, pos + 1)] <- refit$coef
    return(beta.refit)
  }
  else {
    # refit <- lm(y ~ 1)
    beta.refit[1] <- mean(y)
    return(beta.refit)
  }
}


update.beta.new <- function(x, y, beta) {
  pos <- which(abs(beta[-1]) > 1e-10)
  beta.refit <- rep(0, length(beta))
  if (length(pos) >= 1) {
    x_pos <- x[, pos, drop = FALSE]
    # OLS first
    refit_ols <- lm(y ~ x_pos)
    coef_ols <- refit_ols$coef
    if (any(is.na(coef_ols))) {
      # OLS failed，fallback to Ridge
      # cat("Warning: OLS refit produced NA coefficients (singular matrix?). Falling back to Ridge.\n")
      refit_ridge <- glmnet::glmnet(x_pos, y, alpha = 0, lambda = 1e-6, standardize = FALSE)
      beta.refit[1] <- refit_ridge$a0
      beta.refit[pos + 1] <- as.vector(refit_ridge$beta)
    } else {
      # OLS succeed
      beta.refit[c(1, pos + 1)] <- coef_ols
    }
    return(beta.refit)
  } else {
    beta.refit[1] <- mean(y)
    return(beta.refit)
  }
}

########################################
### Calculate mean squared residuals ###
########################################

### Inputs:
# as in update.beta
### Output:
# mean squared residuals

MSR.function <- function(x, y, beta) {
  eta <- cbind(rep(1, nrow(x)), x) %*% beta
  res <- y - eta
  return(sum(res * res)/nrow(x))
}

#######################################################
### pLASSO in linear regression by cross validation ###
#######################################################

### This is the main function for linear regression with cross validation.
### It can find each unweighted penalization estimators: LASSO, p, pLASSO.
### It can find each weighted penalization estimators: LASSO-A, p-A, pLASSO-A.
### For details of the above terms, refer to the paper.

### Inputs:
# x: design matrix
# y: response vector
# eta.seq: the sequence of eta, the parameter tuning the prior information [(6) in the paper]
# for LASSO/p/pLASSO, use eta sequence
# for LASSO-A/p-A/pLASSO-A, use 0
# tau.seq: the sequence of tau, the parameter tuning the L1 penalty weights [(35) in the paper]
# for LASSO/p/pLASSO, use 0
# for LASSO-A/p-A/pLASSO-A, use tau sequence
# beta.prior: either the prior estimator or the initial estimator
# for LASSO/p: use a vector with all 0
# for pLASSO: use prior estimator
# for LASSO-A: use LASSO estimator
# for p-A: use prior estimator
# for pLASSO-A: use pLASSO estimator
# penalty.factor: a vector of indicators of which predictors are penalized (yes/no: 1/0)
# cv.group: the index of cross validation groups
# is.refit: is the estimator refitted (yes/no: 1/0)?
### Output:
# the optimal estimator with cross validation

find.beta <- function(x, y, eta.seq, tau.seq, beta.prior, penalty.factor, cv.group, is.refit)
{
  n.group <- length(unique(cv.group))
  y.prior <- cbind(1, x) %*% beta.prior
  log.liklhd.eta.tau <- matrix(0, length(eta.seq), length(tau.seq))
  betas.eta.tau <- array(0, dim = c(ncol(x) + 1, length(eta.seq), length(tau.seq)))

  for(k in 1 : length(eta.seq))
  {
    eta <- eta.seq[k]
    y.tilde <- (y + eta * y.prior)/(1 + eta)
    for(l in 1 : length(tau.seq))
    {
      tau <- tau.seq[l]
      Winv <- 1 + tau * abs(beta.prior[-1])
      x.Winv <- x * (matrix(1, nrow = nrow(x)) %*% Winv)
      fits <- glmnet(x.Winv, y.tilde, standardize = FALSE, penalty.factor = penalty.factor)

      lambda.seq <- fits$lambda
      log.liklhd <- matrix(0, length(lambda.seq), n.group)
      a0s <- fits$a0
      betas <- fits$beta

      for(j in 1 : n.group)
      {
        sel <- (cv.group != j)
        x.Winv.train <- x.Winv[sel, ]
        y.tilde.train <- y.tilde[sel]
        x.tune <- x[!sel, ]
        y.tune <- y[!sel]
        fits <- glmnet(x.Winv.train, y.tilde.train, standardize = FALSE, penalty.factor = penalty.factor)

        for(i in 1 : length(lambda.seq))
        {
          beta <- coef(fits, s = lambda.seq[i])
          beta.refit <- rep(0, length(beta))
          if(is.refit) {
            temp <- update.beta(x.Winv.train, y.tilde.train, beta)
          }
          else {
            temp <- beta
          }
          beta.refit[1] <- temp[1]
          beta.refit[-1] <- Winv * temp[-1]
          log.liklhd[i, j] <- MSR.function(x.tune, y.tune, beta.refit) * (-1)
          #if(sum(is.na(beta.refit))>0){browser()}
        }
      }

      mean.log.liklhd <- apply(log.liklhd, 1, mean)
      se.log.liklhd <- apply(log.liklhd, 1, sd)/sqrt(n.group)
      print(mean.log.liklhd)
      pos.max <- which.max(mean.log.liklhd)
      # cat("pos.max = ", pos.max, "\n")
      # one standard error rule
      temp <- which(mean.log.liklhd >= mean.log.liklhd[pos.max] - se.log.liklhd[pos.max])
      pos.1se <- temp[1]
      cat("pos.1se = ", pos.1se, "\n")
      log.liklhd.eta.tau[k, l] <- mean.log.liklhd[pos.1se]
      betas.eta.tau[, k, l] <- c(a0s[pos.1se], betas[, pos.1se])
    }
  }

  pos.max.eta.tau <- which(log.liklhd.eta.tau == max(log.liklhd.eta.tau), arr.ind = TRUE)
  if(nrow(pos.max.eta.tau) > 1) {
    cat("pos.max.eta.tau has more than one row", "\n")
    pos.max.eta.tau <- pos.max.eta.tau[nrow(pos.max.eta.tau), ]
  }
  eta <- eta.seq[pos.max.eta.tau[1]]
  tau <- tau.seq[pos.max.eta.tau[2]]

  if(length(eta.seq) > 1 && length(tau.seq) == 1) {
    cat("plasso: selected eta is ", eta, "\n")
  }
  if(length(eta.seq) == 1 && length(tau.seq) > 1) {
    cat("alasso: selected tau is ", tau, "\n")
  }
  if(length(eta.seq) > 1 && length(tau.seq) > 1) {
    cat("plasso + alasso: selected eta is ", eta, " and selected tau is ", tau, "\n")
  }

  y.tilde <- (y + eta * y.prior)/(1 + eta)
  Winv <- 1 + tau * abs(beta.prior[-1])
  x.Winv <- x * (matrix(1, nrow = nrow(x)) %*% Winv)
  beta <- betas.eta.tau[, pos.max.eta.tau[1], pos.max.eta.tau[2]]
  beta.refit <- rep(0, length(beta))
  if(is.refit) {
    temp <- update.beta(x.Winv, y.tilde, beta)
  }
  else {
    temp <- beta
  }
  beta.refit[1] <- temp[1]
  beta.refit[-1] <- Winv * temp[-1]

  return(list(beta = beta.refit,eta = eta, tau = tau))
}

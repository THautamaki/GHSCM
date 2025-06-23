graphical_horseshoe_map <- function(X, fixed_tau = 0, tol = 1e-4, max_iter = 500,
                                    diff_type = "relative", save_lambdas_and_deltas = FALSE,
                                    initial_values = NULL, verbose = 1) {
  start_time <- Sys.time()
  iter <- diff <- 1
  n <- nrow(X)
  p <- ncol(X)
  S <- t(X) %*% X
  diffs <- rep(NA, max_iter)
  if (save_lambdas_and_deltas) {
    lambdas <- deltas <- array(dim = c(p, p, max_iter))
  }
  ind_noi_all <- matrix(0, nrow = p - 1, ncol = p)
  for (i in 1:p) {
    if (i == 1)
      ind_noi <- 2:p
    else if (i == p)
      ind_noi <- 1:(p - 1)
    else
      ind_noi <- c(1:(i - 1), (i + 1):(p))
    ind_noi_all[, i] <- ind_noi
  }
  # Initialise initial values of the algorithm.
  Omega <- Omega_last_iter <- diag(p)
  if (!is.null(initial_values)) {
    if (any(names(initial_values) == "Sigma")) {
      if (matrixcalc::is.positive.definite(initial_values$Sigma) &
          all(dim(initial_values$Sigma) == c(p,p))) {
        Sigma <- initial_values$Sigma
      }
      else {
        cat("Initial value for Sigma is not positive-definite matrix or has wrong dimensions!\nContinue using identity matrix.\n")
        Sigma <- diag(p)
      }
    }
    else {
      cat("Sigma does not found in the initial values! Continue using identity matrix.\n")
      Sigma <- diag(p)
    }
    if (any(names(initial_values) == "Lambda")) {
      if (all(dim(initial_values$Lambda) == c(p,p))) {
        Lambda_sq <- initial_values$Lambda
      }
      else {
        cat("Initial value for Lambda has wrong dimensions!\nContinue using one matrix.\n")
        Lambda_sq <- matrix(1, p, p)
      }
    }
    else {
      cat("Lambda does not found in the initial values! Continue using one matrix.\n")
      Lambda_sq <- matrix(1, p, p)
    }
    if (any(names(initial_values) == "Delta")) {
      if (all(dim(initial_values$Delta) == c(p,p))) {
        Delta <- initial_values$Delta
      }
      else {
        cat("Initial value for Delta has wrong dimensions!\nContinue using one matrix.\n")
        Delta <- matrix(1, p, p)
      }
    }
    else {
      cat("Delta does not found in the initial values! Continue using one matrix.\n")
      Delta <- matrix(1, p, p)
    }
  }
  else {
    Sigma <- diag(p)
    Lambda_sq <- Delta <- matrix(1, p, p)
  }
  tau_sq <- fixed_tau
  while (diff > tol & iter <= max_iter) {
    for (i in 1:p) {
      # Extract matrix blocks into own variables.
      ind_noi <- ind_noi_all[, i]
      Sigma_11 <- Sigma[ind_noi, ind_noi]
      Sigma_12 <- Sigma[ind_noi, i]
      Sigma_22 <- Sigma[i, i]
      lambda_sq_12 <- Lambda_sq[ind_noi, i]
      delta_12 <- Delta[ind_noi, i]
      s_12 <- S[ind_noi, i]
      s_22 <- S[i, i]
      # Calculate gamma, omega_12, omega_22, lambda_12^2 and delta_12.
      gamma <- (n / 2) / (s_22 / 2)
      Omega_11_inv <- Sigma_11 - Sigma_12 %*% t(Sigma_12) / Sigma_22
      C_inv <- s_22 * Omega_11_inv + diag(1 / (lambda_sq_12 * tau_sq))
      beta <- -solve(C_inv, s_12)                  # Same as -C %*% s_12
      omega_12 <- beta
      omega_22 <- gamma + t(beta) %*% Omega_11_inv %*% beta
      lambda_scale <- 1 / delta_12 + omega_12^2 / (2 * tau_sq)
      lambda_sq_12 <- lambda_scale / 2
      delta_scale <- 1 + 1 / lambda_sq_12
      delta_12 <- delta_scale / 2
      # Update Omega, Sigma, Lambda^2 and Delta matrices.
      Omega[i, i] <- omega_22
      Omega[i, ind_noi] <- omega_12
      Omega[ind_noi, i] <- omega_12
      Omega_inv_temp <- Omega_11_inv %*% beta
      Sigma[ind_noi, ind_noi] <- Omega_11_inv + Omega_inv_temp %*% t(Omega_inv_temp) / gamma
      Sigma_12 <- -Omega_inv_temp / gamma
      Sigma[ind_noi, i] <- Sigma_12
      Sigma[i, ind_noi] <- Sigma_12
      Sigma[i, i] <- 1 / gamma 
      Lambda_sq[i, ind_noi] <- lambda_sq_12
      Lambda_sq[ind_noi, i] <- lambda_sq_12
      Delta[i, ind_noi] <- delta_12
      Delta[ind_noi, i] <- delta_12
    }
    if (diff_type == "relative") {
      diff <- norm(Omega - Omega_last_iter, type = "F") / norm(Omega_last_iter, type = "F")
    }
    else {
      diff <- norm(Omega - Omega_last_iter, type = "F")
    }
    diffs[iter] <- diff
    if (verbose > 1) {
      lap_time <- Sys.time()
      elap_time <- as.numeric(difftime(lap_time, start_time, unit = "s"))
      cat("Iteration: ", iter, ". Elapsed time: ", round(elap_time, 2), " s. Difference: ",
          round(diff, 6), "\n", sep = "")
    }
    else if (verbose > 0 & iter %% 10 == 0) {
      lap_time <- Sys.time()
      elap_time <- as.numeric(difftime(lap_time, start_time, unit = "s"))
      cat("Iteration: ", iter, ". Elapsed time: ", round(elap_time, 2), " s. Difference: ",
          round(diff, 6), "\n", sep = "")
    }
    Omega_last_iter <- Omega
    if (save_lambdas_and_deltas) {
      lambdas[,,iter] <- Lambda_sq
      deltas[,,iter] <- Delta
    }
    iter <- iter + 1
  }
  if (save_lambdas_and_deltas) {
    lambdas <- lambdas[,, 1:(iter-1)]
    deltas <- deltas[,, 1:(iter-1)]
  }
  diffs <- diffs[(1:iter-1)]
  
  kappa <- 1 / (1 + n * tau_sq * Lambda_sq)
  theta <- 1 - kappa
  theta[theta >= 0.5] <- 1
  theta[theta < 0.5] <- 0
  diag(theta) <- 0
  if (verbose >= 0) {
    end_time <- Sys.time()
    elap_time <- as.numeric(difftime(end_time, start_time, unit = "s"))
    cat("Total iterations: ", iter - 1, ". Elapsed time: ", round(elap_time, 2),
        " s. Final difference: ", diff, "\n", sep = "")
  }
  if (save_lambdas_and_deltas) {
    return(list(Sigma_est = Sigma, Omega_est = Omega, Theta_est = theta, Kappa = kappa, Lambdas = lambdas,
                tau = tau_sq, Deltas = deltas, diffs = diffs, iters = iter - 1, tot_time = elap_time))
  }
  else {
    return(list(Sigma_est = Sigma, Omega_est = Omega, Theta_est = theta, Kappa = kappa,
                tau = tau_sq, diffs = diffs, iters = iter - 1, tot_time = elap_time))
  }
}

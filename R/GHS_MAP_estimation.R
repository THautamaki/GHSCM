
GHS_MAP_estimation <- function(X, p0 = 0, tolerance = 1e-4, max_iterations = 500, use_Cpp = TRUE,
                               save_lambdas_and_deltas = FALSE, initial_values = NULL, verbose = 0) {
  n <- nrow(X)
  p <- ncol(X)
  if (p0 == 0) {
    p0 <- p - 1
  }
  if (p < 200) {
    tau_f <- (p0/(p*(p-1)/2)) * (50 * sqrt(p0) / (n * sqrt(n)))
  }
  else {
    tau_f <- (p0/(p*(p-1)/2)) * ((p/4) * sqrt(p0) / (n * sqrt(n)))
  }
  if (use_Cpp) {
    GHSGEM_est <- graphical_horseshoe_map_Cpp(X, tol = tolerance, max_iter = max_iterations,
                                              fixed_tau = tau_f, verbose = verbose)
    GHSGEM_est$diffs <- GHSGEM_est$diffs[1:GHSGEM_est$iters]
  }
  else {
    GHSGEM_est <- graphical_horseshoe_map(X, tol = tolerance, max_iter = max_iterations,
                                          fixed_tau = tau_f, initial_values = initial_values,
                                          save_lambdas_and_deltas = save_lambdas_and_deltas,
                                          verbose = verbose)
  }
  return(GHSGEM_est)
}
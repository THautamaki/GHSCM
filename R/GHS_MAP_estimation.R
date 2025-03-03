
#' Calculate MAP estimates for the precision and covariance matrices and for the network.
#' 
#' This function calculates MAP estimates for the precision matrix and network using generalised
#' expectation-maximisation (GEM) algorithm with graphical horseshoe (GHS) prior. The input data
#' should follow a multivariate normal distribution with Var(\eqn{X_i}) = 1, for \eqn{\forall \, i}.
#'
#' @param X The data matrix with dimensions \code{n} by \code{p}, where \code{n} = sample size and
#'   \code{p} = number of variables.
#' @param p0 A prior assumption of number of connections in the network. The default values is 0,
#'   which means that p - 1 is used.
#' @param tolerance Convergence tolerance of the algorithm. The default value is 1e-4.
#' @param max_iterations Maximum iterations of the algorithm, if convergence tolerance is not
#'   reached. The default value is 500 iterations.
#' @param use_Cpp Boolean, which used to select C++ or R version of algorithm. The default value is
#'   \code{TRUE} (C++ version).
#' @param save_lambdas_and_deltas Boolean, which used to select if values of the Lambda and Delta
#'   matrices are saved. The default value is \code{FALSE}. Only used with the R version of
#'   algorithm (\code{use_Cpp = FALSE}).
#' @param initial_values Named list of the initial values of the algorithm. Must contain names
#'   \code{Sigma}, \code{Lambda} and/or \code{Delta}. If wrong values used, reverts back to default
#'   initial values. Only used with the R version of algorithm (\code{use_Cpp = FALSE}).
#' @param verbose Numeric with possible values of 0, 1 and 2. 0 (default) prints only information of
#'   final iteration. 1 prints information every 10th iteration, and 2 prints information every
#'   iteration.
#'
#' @return
#' A list, which contains following objects:
#' \item{Sigma_est}{
#'    The estimate for the \code{p} by \code{p} covariance matrix.
#' }
#' \item{Omega_est}{
#'    The estimate for the \code{p} by \code{p} precision matrix.
#' }
#' \item{Theta_est}{
#'    The \code{p} by \code{p} adjacency matrix of the network estimate.
#' }
#' \item{Kappa}{
#'    The \code{p} by \code{p} Kappa matrix, which was used to construct sparse network estimate.
#' }
#' \item{tau}{
#'    The actual value of the tau^2 parameter.
#' }
#' \item{diffs}{
#'    The convergence differences of each iteration.
#'    }
#' \item{iters}{
#'    The number of algorithm iterations.
#' }
#' \item{tot_time}{
#'    The total runtime of the algorithm in seconds.
#' }
#' If \code{save_lambdas_and_deltas = TRUE} and the R version of algorithm is in use (\code{use_Cpp = FALSE},
#' then return also:
#' \item{Lambdas}{
#'    \code{p} by \code{p} by \code{iters} array, which contains the Lambda matrices of
#'    each iteration.
#' }
#' \item{Deltas}{
#'    \code{p} by \code{p} by \code{iters} array, which contains the Delta matrices of
#'    each iteration.
#' }
#'   
#' @export
#'
#' @examples
#' sim <- huge::huge.generator(n = 200, d = 100, graph = "scale-free") # Generate data
#' map <- GHS_MAP_estimation(sim$data, verbose = 1)                    # Run algorithm
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
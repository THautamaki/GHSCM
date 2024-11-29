
GHS_MAP_estimate <- function(X, p0 = 0, verbose = 0, tol = 1e-4, max_iterations = 200) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (p0 == 0) {
    p0 <- p - 1
  }
  tau_f <- (p0/n)^1.5*(50/(p*(p-1)/2))
  map <- graphical_horsesoe_map_Cpp(X, verbose = verbose, tol = tol, max_iter = max_iterations)
  return(map)
}
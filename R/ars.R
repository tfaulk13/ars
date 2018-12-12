#' Adaptive Rejection Sampler
#'
#' Sample points from a distribution using an adaptive rejection sampler.
#'
#' @param n the number of samples that you'd like taken from the distribution.
#' @param f the function used to determine the distribution from which you'll
#'   sample
#' @param dfunc the derivative of the function.
#' @param D_min the minimum boundary limit of the distribution.
#' @param D_max the maximum boundary limit of the distribution.
#' @param verbose Whether you want R to print out the process as it goes through
#'   ARS. This can be useful for longer computations.
#'
#' @details Adaptive rejection sampling can be beneficial when sampling from
#' certain distributions by reducing the number of evaluations that must occur.
#' It reduces the evaluations by assuming log-concavity for any input function,
#' and because it doesn't need to update for new envelope and squeeze functions
#' every iteration. ARS is particularly useful in Gibbs Sampling, where
#' calculations are complicated but usually log-concave.
#'
#'
#' @examples
#' # testing on normal distribution
#' res_samples <- ars(n = 1000, f = dnorm)
#'
#'
#' # Check variance and mean
#' mean(res_samples)
#' var(res_samples)
#'
ars <- function(n, f, starting_points, dfunc, interval) {

  if (missing(n))
    stop("Please provide the number of samples 'n' desired.")
  if (!is.function(f))
    stop("Please provide an R function as input for 'f'.")

  log_f <- function(x) log(f(x))

  if (missing(dfunc)) {
    dfunc <- function(x) fderiv(f, x)/(f(x) + .Machine$double.eps)
  }

  # Restrict the density to a smaller interval containing the overall information of the density.
  if (missing(interval)) {
    support_limits <- set_support_limit(f)
    D_min <- support_limits$D_min
    D_max <- support_limits$D_max
  } else {
    assertthat::are_equal(length(interval), 2)
    testthat::expect_type(interval, "double")
    D_min <- interval[1]
    D_max <- interval[2]
    if (log_f(D_min) == -Inf || log_f(D_max) == -Inf) {
      stop("Please provide a more compact interval.")
    }
  }





  # Restrict the density to a smaller interval containing the overall information of the density.
  # if (missing(interval)) {
  #   support_limits <- set_support_limit(log_f, dfunc, D_min = -1E8, D_max = 1E8)
  #   D_min <- support_limits$D_min
  #   D_max <- support_limits$D_max
  # } else {
  #   assertthat::are_equal(length(interval), 2)
  #   testthat::expect_type(interval, "double")
  #   D_min <- interval[1]
  #   D_max <- interval[2]
  #   if (log_f(D_min) == -Inf || log_f(D_max) == -Inf) {
  #     stop("Pleasee provide a more compact interval.")
  #   }
  # }


  iteration <- 0
  samples <- c()

  # Come up with starting points if not given.
  if (missing(starting_points)) {
    T_current <- find_init_points(f = log_f, dfunc = dfunc, D_min = D_min, D_max = D_max)
  } else {
    if (!is.numeric(starting_points))
      stop("The 'starting_points' input should be a vector of abscissae.")
    if (length(starting_points) < 2)
      stop("The 'starting_points' input should contain at least 2 elements.")
    T_current <- starting_points
  }

  # Find Tangent intersections points.
  z_points <- sort(c(D_min, find_tangent_intercept(log_f, dfunc, T_current), D_max))

  # Define the upper and lower envelopes, and the sampling function.
  helper_funcs <- get_updated_functions(T_ = T_current, z_pts = z_points, func = log_f, dfunc = dfunc, D_min = D_min, D_max = D_max)
  u_current <- helper_funcs$upper
  s_current <- helper_funcs$s
  l_current <- helper_funcs$lower
  integrals_s <- helper_funcs$s_integrals

  while (length(samples) <= n) {
    x_star <- sample_x_star(s_function = s_current, s_integrals = integrals_s, z_pts = z_points)

    w <- runif(1)

    ## 1st squeezing test.
    diff_upper_lower <- u_current(x_star) - l_current(x_star)
    if (diff_upper_lower < -1E-5) stop("Function is not log concave. Adaptive rejection sampling won't give a proper sample for this distribution.")
    if ( w <= exp(-diff_upper_lower) ) {
      samples <- c(samples, x_star)
    } else {
      diff_upper_f <- u_current(x_star) - log_f(x_star)
      if(diff_upper_f < -1E-5) stop("Function is not log concave. Adaptive rejection sampling won't give a proper sample for this distribution.")
      if (w <= exp(-diff_upper_f)) {
        samples <- c(samples, x_star)
      }

      # Update step.
      T_current <- sort(c(T_current, x_star))
      iteration <- iteration + 1

      z_points <- (sort(c(D_min, find_tangent_intercept(log_f, dfunc, T_current), D_max)))
      helper_funcs <- get_updated_functions(T_ = T_current, z_pts = z_points, func = log_f, dfunc = dfunc, D_min = D_min, D_max = D_max)
      u_current <- helper_funcs$upper
      s_current <- helper_funcs$s
      l_current <- helper_funcs$lower
      integrals_s <- helper_funcs$s_integrals
    }
  }

  return(samples)
}

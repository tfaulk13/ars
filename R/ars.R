#' Adaptive Rejection Sampler
#'
#' Sample points from a distribution using an adaptive rejection sampler. Based upon:
#' Gilks, W., & Wild, P. (1992). Adaptive Rejection Sampling for Gibbs Sampling. Journal of the Royal Statistical Society.
#'
#' @param n the number of samples that you'd like taken from the distribution.
#' @param f the function used to determine the distribution from which you'll
#'   sample
#' @param dfunc the derivative of the function.
#' @param interval the boundary limit of the distribution.
#' @param starting_points the starting points used for the initialization step of the ars algorithm.
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
#' # specifying a domain interval
#' res_samples <- ars(n, function(x) dnorm(x, mean = 3, sd = 5), interval = c(-20, 25))
#'
#'# specifying starting points (check the paper in references for more details)
#' res_samples <- ars(n, dunif, interval = c(0,1), starting_points = c(0.2, 0.8))
#'
#' @references
#' Gilks, W., & Wild, P. (1992). Adaptive Rejection Sampling for Gibbs Sampling. Journal of the Royal Statistical Society.
#'
#' @export
#'
ars <- function(n, f, starting_points, dfunc, interval) {

  if (missing(n))
    stop("Please provide the number of samples 'n' desired.")
  if (!is.function(f))
    stop("Please provide an R function as input for 'f'.")

  log_f <- function(x) log(f(x))

  if (missing(dfunc)) {
    dfunc <- function(x) pracma::fderiv(f, x) / (f(x) + .Machine$double.eps)
  }

  # Restrict the density to a smaller interval containing the overall information of the density.
  if (missing(interval)) {
    support_limits <- get_support_limit(f)
    D_min <- support_limits$D_min
    D_max <- support_limits$D_max
  } else {
    assertthat::are_equal(length(interval), 2)
    testthat::expect_type(interval, "double")
    testthat::expect_equal(length(interval), 2)
    D_min <- interval[1] + .Machine$double.eps^0.25
    D_max <- interval[2] - .Machine$double.eps^0.25
    if (log_f(D_min) == -Inf || log_f(D_max) == -Inf) {
      stop("Please provide a more compact interval, where the density does not evaluate to 0.")
    }
  }

  samples <- c()

  # Come up with starting points if not given.
  if (missing(starting_points)) {
    T_current <- find_init_points(f = log_f, dfunc = dfunc, D_min = D_min, D_max = D_max)
  } else {
    if (!is.numeric(starting_points))
      stop("The 'starting_points' input should be a vector of abscissae.")
    if (length(starting_points) < 2)
      stop("The 'starting_points' input should contain at least 2 elements.")
    if (Inf %in% abs(log_f(starting_points)))
      stop("Please give starting points where the density can be evaluated with a minimum precision.")
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

  while (length(samples) < n) {
    # sample an element from s_current.
    x_star <- sample_x_star(s_function = s_current, s_integrals = integrals_s, z_pts = z_points)

    w <- runif(1)

    ## 1st squeezing test.
    diff_upper_lower <- u_current(x_star) - l_current(x_star)
    if (diff_upper_lower < -1E-5) stop("Function is not log concave. Adaptive rejection sampling won't give a proper sample for this distribution.")
    if ( w <= exp(-diff_upper_lower) ) {
      samples <- c(samples, x_star)
    } else {
      # 2nd squeezing test
      diff_upper_f <- u_current(x_star) - log_f(x_star)
      if(diff_upper_f < -1E-5) stop("Function is not log concave. Adaptive rejection sampling won't give a proper sample for this distribution.")
      if (w <= exp(-diff_upper_f)) {
        samples <- c(samples, x_star)
      }

      # Update step.
      T_current <- sort(c(T_current, x_star))

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

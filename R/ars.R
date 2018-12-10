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
#' h = function(y) {
#' return(-log(2 * pi * sigma ^ 2) / 2 - (y - mu) ^ 2 / (2 * sigma ^ 2))
#' }
#'
#' # derivative of h
#' dh = function(y){
#'   return(- (y - mu)/sigma^2)
#' }
#'
#' res_samples <- ars(n = 1000, f = h, dfunc = dh, D_min = -10, D_max = 30)
#'
#'
#' # Check variance and mean
#' mean(res_samples)
#' var(res_samples)
#'
ars <- function(n, f, dfunc, D_min, D_max, verbose = FALSE) {
  ## Need to work on this. When D_min or D_max = +- Inf problems occur.
  if (D_min == -Inf) {
    D_min = -10^8
  }
  if (D_max == Inf) {
    D_max = 10^8
  }

  iteration <- 0
  samples <- c()
  T_current <- find_init_points(f, dfunc, D_min, D_max)

  # tangent intersections points.
  z_points <- sort(c(D_min, find_tangent_intercept(f, dfunc, T_current), D_max))

  helper_funcs <- get_updated_functions(T_ = T_current, z_pts = z_points, func = f, dfunc = dfunc, D_min = D_min, D_max = D_max)
  u_current <- helper_funcs$upper
  s_current <- helper_funcs$s
  l_current <- helper_funcs$lower


  while (length(samples) <= n) {
    # cat("T: ", T_current, "\n")
    # cat("samples: ", samples, "\n")
    integrals_s <- sapply(X = seq(1, length(z_points) -1), FUN = function(i) integrate(s_current, lower = z_points[i], upper = z_points[i + 1])$value)
    idx <- sample(seq(1, length(z_points) -1), size = 1, prob = integrals_s)

    unnormalized_cdf_s <- function(x) {
      if (x <= z_points[idx]) {
        return(0)
      } else if (x >= z_points[idx + 1]) {
        return(integrals_s[idx])
      } else {
        integrate(s_current, lower = z_points[idx], upper = x)$value
      }
    }

    quantile <- runif(1, min = 0, max = integrals_s[idx])
    x_star <- uniroot(f = function(x) unnormalized_cdf_s(x) - quantile, interval = c(z_points[idx], z_points[idx + 1]))$root # sample from s_current

    assertthat::are_equal(x_star <= z_points[idx + 1] && x_star >= z_points[idx], TRUE)
    w <- runif(1)

    ## 1st squeezing test.
    if ( w <= exp(l_current(x_star) - u_current(x_star)) ) {
      samples <- c(samples, x_star)
      # cat("x_star accepted at 1st test. Samples: ", samples, "\n")
    } else {
      if (w <= exp(f(x_star) - u_current(x_star))) {
        # if (verbose == TRUE) {
        #   cat("The point ", x_star, " has been added to the abscissae points at iteration ", iteration, "\n")
        # }
        samples <- c(samples, x_star)
        #cat("x_star accepted at 2nd test. Samples: ", samples, "\n")
      }

      # Update step
      T_current <- sort(c(T_current, x_star))
      iteration <- iteration + 1

      z_points <- sort(c(D_min, find_tangent_intercept(f, dfunc, T_current), D_max))
      helper_funcs <- get_updated_functions(T_ = T_current, z_pts = z_points, func = f, dfunc = dfunc, D_min = D_min, D_max = D_max)
      u_current <- helper_funcs$upper
      s_current <- helper_funcs$s
      l_current <- helper_funcs$lower
    }
  }

  return(samples)
}

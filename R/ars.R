## Loading libraries ----------------------------------------------------------
library(devtools)
library(roxygen2)

## find_tangent_intercept -----------------------------------------------------
#
# This function finds the intercept of the tangent lines at the beginning
#
# Inputs:
# f = a function provided by the user.
# df = the derivative of the function.
# x_pts = points selected initially to find the tangent line.
#
# Outputs:
# z = the intercepts of the tangent lines
#
find_tangent_intercept = function(f, df, x_pts){
  x0 = head(x_pts, n = -1)
  x1 = tail(x_pts, n = -1)
  z = (f(x1) - f(x0) - x1 * df(x1) + x0 * df(x0)) / (df(x0) - df(x1))

  assertthat::are_equal(length(x_pts), length(z) + 1)

  return(z)
}

## get_upper ------------------------------------------------------------------
#
# This function calculates the upper bound
#
# Inputs:
# T_ = the abcissa of the function.
# z_pts = the intercepts of the tangent line.
# func = a function provided by the user.
# dfunc = the derivative of the function.
#
# Output:
# the upper bound.
#
get_upper <- function(T_, z_pts, func, dfunc) {
  return(
    function(x) {
      #cat(z_pts, "\n \n")
      z_idx <- findInterval(x = x, vec = z_pts)
      return(func(T_[z_idx]) + (x - T_[z_idx]) * dfunc(T_[z_idx]))
    }
  )

}

## get_lower ------------------------------------------------------------------
#
#
#
# Inputs:
# T_ = the abcissa of the function.
# func = a function provided by the user.
# dfunc = the derivative of the function.
#
# Output:
# the lower bound.
#
get_lower <- function(T_, func, dfunc) {
  return(
    function(x) {
      res <- c()
      x_idx <- findInterval(x, T_)
      for (i in seq(1, length(x))) {
        if (x[i] %in% T_) {
          res[i] <- func(T_[x_idx[i]])
        } else if (x_idx[i] == 0 || x_idx[i] == length(T_)) {
          res[i] <- -Inf
        } else {
          numerator <- (T_[x_idx[i] + 1] - x[i]) * func(T_[x_idx[i]]) + (x[i] - T_[x_idx[i]]) * func(T_[x_idx[i] + 1])
          denominator <- T_[x_idx[i] + 1] - T_[x_idx[i]]
          res[i] <- numerator / denominator
        }
      }
      return(res)

      # numerator <- (T_current[x_idx + 2] - x)*f_(T_current[x_idx + 1]) + (x - T_current[x_idx + 1])*f_(T_current[x_idx + 2])
      # denominator <- T_current[x_idx + 2] - T_current[x_idx + 1]
      # return( numerator/denominator)
    }
  )
}

## Get updated functions ------------------------------------------------------
#
# This function updates the functions as neded when ARS needs new inputs.
#
# Inputs:
# T_ = the abcissa of the function.
# z_pts = the intercepts of the tangent line.
# func = a function provided by the user.
# dfunc = the derivative of the function.
# D_min = the minimum boundary limit.
# D_max = the maximimum boundary limit of the distribution.
#
# Outputs:
# An updated upper bound, lower bound, and new starting points for the next
# iteration of the function.
#
get_updated_functions <- function(T_, z_pts, func, dfunc, D_min, D_max) {

  u_current <- get_upper(z_pts = z_pts, T_ = T_, func = func, dfunc = dfunc)

  s_current <- function(x) {
    exp(u_current(x)) / integrate(function(x) exp(u_current(x)), lower = D_min, upper = D_max)$value
  }
  l_current <- get_lower(T_ = T_, func = func, dfunc = dfunc)

  return(list("upper" = u_current,
              "lower" = l_current,
              "s" = s_current
  )
  )
}




#' Adaptive Rejection Sampler
#'
#' Sample points from a distribution using an adaptive rejection sampler.
#'
#' @param n the number of samples that you'd like taken from the distribution.
#' @param f the function used to determine the distribution from which you'll
#'   sample
#' @param dfunc the derivative of the function.
#' @param T_start The starting points to begin ARS.
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
#' res_samples <- ars(n = 1000, f = h, dfunc = dh, T_start = c(-3, 0, 3), D_min
#' = -10, D_max = 30)
#'
#'
#' # Check variance and mean
#' mean(res_samples)
#' var(res_samples)
#'
ars <- function(n, f, dfunc, T_start, D_min, D_max, verbose = FALSE) {
  ## Need to work on this. When D_min or D_max = +- Inf problems occur.
  f_ <- function(x) {
    res <- c()
    for (i in seq(1, length(x))) {
      if (x[i] == Inf || x[i] == -Inf) {
        res[i] <- -Inf
      } else {
        res[i] <- f(x[i])
      }
    }
    return(res)
  }

  df_ <- function(x) {
    res <- c()
    for (i in seq(1, length(x))) {
      if (x[i] == Inf ) {
        res[i] <- dfunc(.Machine$double.xmax)
      } else if ( x[i] == -Inf) {
        res[i] <- dfunc(.Machine$double.xmin)
      } else {
        res[i] <- dfunc(x[i])
      }
    }
    return(res)
  }

  iteration <- 0
  samples <- c()
  T_current <- T_start

  # tangent intersections points.
  z_points <- sort(c(D_min, find_tangent_intercept(f, dfunc, T_current), D_max))

  helper_funcs <- get_updated_functions(T_ = T_current, z_pts = z_points, func = f_, dfunc = df_, D_min = D_min, D_max = D_max)
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
      helper_funcs <- get_updated_functions(T_ = T_current, z_pts = z_points, func = f_, dfunc = df_, D_min = D_min, D_max = D_max)
      u_current <- helper_funcs$upper
      s_current <- helper_funcs$s
      l_current <- helper_funcs$lower
    }
  }

  return(samples)
}





## Test with h

params.r = 2
params.m = 10
params.mu = 0
params.sig2 = 1

h = function(y){
  v = params.r*y - params.m * log(1+exp(y)) - (y-params.mu)^2/(2*params.sig2) # plus normalizing const
  return(v)
}
## derivative of h
dh = function(y){
  params.r - params.m * exp(y) / (1 + exp(y)) - (y-params.mu)/params.sig2
}

res_samples <- ars(n = 10, f = h, dfunc = dh, T_start = c(-3, 0, 3), D_min = -20, D_max = 20)


## Test with normal distribution
mu = 10
sigma = 2

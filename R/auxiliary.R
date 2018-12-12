## find_init_points -----------------------------------------------------------
#
# This function finds the initial starting points for ARS. Used when the user does give it as an input
#
# Inputs:
# f = a function provided by the user.
# dfunc = the derivative of the function.
# D_min = the minimum boundary limit.
# D_max = the maximimum boundary limit of the distribution.
#
# Outputs:
# The starting points for ARS
#
find_init_points <- function(f, dfunc, D_min, D_max) {
  # case that the func is monotonic increasing/decreasing
  if (sign(dfunc(D_min)) == sign(dfunc(D_max)) ) {
    X_init <- c(D_min, (D_min + D_max) / 2, D_max)
  } else{# the case that there is a max in the range
    gap <- (D_max - D_min) / 200
    max <- optimize(f = f, interval = c(D_min, D_max), lower = D_min, upper = D_max, maximum = TRUE)$maximum
    right_pt <- (max + gap)
    left_pt <- (max - gap)
    X_init <- c(left_pt, max, right_pt)
  }
  return(sort(X_init))
}

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
# the upper bound function.
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
# the lower bound function.
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
    }
  )
}

## get_updated_functions ------------------------------------------------------
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

  integral_u <- integrate(function(x) exp(u_current(x)), lower = D_min, upper = D_max)$value
  s_current <- function(x) {
    res <- exp( u_current(x) - log(integral_u))
    #res[which(is.na(res))] <- 0
    return(res)
  }
  l_current <- get_lower(T_ = T_, func = func, dfunc = dfunc)

  s_integrals <- sapply(X = seq(1, length(z_pts) -1), FUN = function(i) integrate(s_current, lower = z_pts[i], upper = z_pts[i + 1])$value)

  return(list("upper" = u_current,
              "lower" = l_current,
              "s" = s_current,
              "s_integrals" = s_integrals)
  )
}

## sample_x_star ------------------------------------------------------
#
# This function samples an element from the piecewise exponential density distribution s_function.
#
# Inputs:
# s_function = the density to sample from.
# z_pts = the piecewise subdivisions of s_function.
# s_integrals = integral of s_function over its subdivisions.
# Outputs:
# A randomly generated value from the s_function density distribution.
#
sample_x_star <- function(s_function, z_pts, s_integrals) {
  if (missing(s_integrals)) {
    # Applying the closed form solution of the integral instead, computed from the Adaptive Rejection Sampling paper would be quicker, but less explicit.
    s_integrals <- sapply(X = seq(1, length(z_pts) -1), FUN = function(i) integrate(s_function, lower = z_pts[i], upper = z_pts[i + 1])$value)
  }

  # Pick the interval of the piecewise exponential distribution we are sampling from.
  idx <- sample(seq(1, length(z_pts) -1), size = 1, prob = s_integrals)

  unnormalized_cdf_s <- function(x) {
    if (x <= z_pts[idx]) {
      return(0)
    } else if (x >= z_pts[idx + 1]) {
      return(s_integrals[idx])
    } else {
      integrate(s_function, lower = z_pts[idx], upper = x)$value
    }
  }

  quantile <- runif(1, min = 0, max = s_integrals[idx])
  x_star <- uniroot(f = function(x) unnormalized_cdf_s(x) - quantile, interval = c(z_pts[idx], z_pts[idx + 1]))$root # sample from s_current

  assertthat::are_equal(x_star <= z_pts[idx + 1] && x_star >= z_pts[idx], TRUE)
  return(x_star)

}


## get_support_limit ------------------------------------------------------
#
# This function defines a domain for the density function so that the main mass of the distribution is within the domain bounds.
#
# Inputs:
# f = the density we are trying to define a domain for.
# Outputs:
# A list containing two elements, the minimum of the domain, and its maximum.
get_support_limit <- function (f) {
  min <- -1E3
  max <- 1E3
  lower_quantile <- log(1E-6)
  upper_quantile <- log(1 - 1E-6)
  cdf <- function(x) {
    norm <- integrate(f, lower = min, upper = max)$value
    if (norm == 0) {
      stop("The density integrates to 0 because of the 'integrate' function behavior on a large interval for a condensed function.
             Give a small interval as input to the ars function and try again.")
    }
    res <- vector(length = length(x))
    for (i in seq(1, length(x))) {
      res[i] <- log(integrate(f, lower = min, upper = x[i])$value) - log(norm)
    }
    return(res)
  }
  safety <- 10
  D_min <- rootSolve::uniroot.all(f = function(x) cdf(x) - lower_quantile, lower = min, upper =  max)[1] - safety
  D_max <- rootSolve::uniroot.all(f = function(x) cdf(x) - upper_quantile, lower = min, upper =  max)[1] + safety

  return(list("D_min" = ifelse(is.na(D_min), -100, D_min),
              "D_max" = ifelse(is.na(D_max), 100, D_max)))
}

## find_init_points -----------------------------------------------------------
#
# This function finds the initial starting points for ARS
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
# s_integrals =
# Outputs:
# A randomly generated valuee from the s_function density distribution.
#
sample_x_star <- function(s_function, z_pts, s_integrals) {
  if (missing(s_integrals)) {
    # Apply the closed form solution of the integral instead ?
    s_integrals <- sapply(X = seq(1, length(z_pts) -1), FUN = function(i) integrate(s_function, lower = z_pts[i], upper = z_pts[i + 1])$value)
  }
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


## set_support_limit ------------------------------------------------------
#
# This function samples an element from the piecewise exponential density distribution s_function.
#
# Inputs:
# s_function = the density to sample from.
# z_pts = the piecewise subdivisions of s_function.
# s_integrals =
# Outputs:
# A randomly generated valuee from the s_function density distribution.
#
# set_support_limit <- function (f, dfunc, D_min = -Inf, D_max = Inf) {
#   clip_value <- -20
#   eps <- 1E-3
#
#   while (f(D_min) < clip_value) {
#     # Newton method adaptation.
#     D_min <- D_min - (f(D_min) - clip_value)/((dfunc(D_min)) + eps)
#   }
#
#   while (f(D_max) < clip_value) {
#     D_max <- D_max + (max(f(D_max), -1E6) - clip_value)/((dfunc(D_max)) + eps)
#   }
#
#   # D_min <- -1E8
#   # while (f(D_min) < clip_value) {
#   #   D_min <- D_min - (f(D_min) - clip_value)/(abs(dfunc(D_min)) + eps)
#   #
#
#   # roots <- rootSolve::uniroot.alluniroot(function(x) f(x) - clip_value, lower = D_min, upper = max)
#   # D_max <- uniroot(function(x) f(x) - clip_value, lower = max, upper = D_max)
#
#   # # We consider that if f(x) is inferior than clip_value,
#   # # then there's no chance of sampling x, as the density at x would correspond to exp(-20).
#   # while (f(D_min) - clip_value > 0) {
#   #   # Newton method adaptation.
#   #   D_min_update <- D_min - (f(D_min) - clip_value)/(abs(dfunc(D_min)) + eps)
#   #   if (f(D_min_update) == -Inf) {
#   #     break
#   #   }
#   #   D_min <- D_min_update
#   # }
#   #
#   #
#   # D_max <- 0
#   # while (f(D_max) > clip_value) {
#   #   D_max_update <- D_max + (f(D_max) - clip_value)/(abs(dfunc(D_max)) + eps)
#   #   if (f(D_max_update) == -Inf) {
#   #     break
#   #   }
#   #   D_max <- D_max_update
#   # }
#   #
#
#   return(list("D_min" = D_min,
#               "D_max" = D_max))
# }

set_support_limit <- function (f) {
  min <- -1E3
  max <- 1E3
  lower_quantile <- 1E-6
  upper_quantile <- 1 - 1E-6
  cdf <- function(x) {
    norm <- integrate(f, lower = min, upper = max)$value
    res <- vector(length = length(x))
    for (i in seq(1, length(x))) {
      res[i] <- integrate(f, lower = min, upper = x[i])$value
    }
    return(res)
  }
  safety <- 10
  D_min <- rootSolve::uniroot.all(f = function(x) cdf(x) - lower_quantile, lower = min, upper =  max)[1] - safety
  D_max <- rootSolve::uniroot.all(f = function(x) cdf(x) - upper_quantile, lower = min, upper =  max)[1] + safety

  return(list("D_min" = D_min,
              "D_max" = D_max))
}

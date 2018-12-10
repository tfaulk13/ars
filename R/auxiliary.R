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
find_init_points <- function(f, dfunc, D_min, D_max){
  # case that the func is monotonic increasing/decreasing
  if ((dfunc(f, D_max) * dfunc(f, D_min)) > 0){
    X_init <- c(D_min, (D_min + D_max) / 2, D_max)
  } else{# the case that there is a max in the range
    gap <- (D_max - D_min) / 200
    max <- optimize(f = h, interval = c(D_min, D_max), lower = D_min, upper = D_max, maximum = TRUE)$maximum
    right_pt <- (max + gap)
    left_pt <- (max - gap)
    X_init <- c(left_pt,max,right_pt)
  }
  return(X_init)
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

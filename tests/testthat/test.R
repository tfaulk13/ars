# ARS Function testing
# "Test-Time" Tests

# Testing 4-6 distributions (quadratic, logistic, exponential, normal)
# Trying different ranges (-inf to inf)

## Log Concavity
context("Tests for Log Concavity")

# Test 1: quadratic function
test_that("tangent lines intercept at zero for quadratic function", {
  f <- x ^ 2
  df <- 2 * x
  x_pts <- c(-1,1) # what value(s) go in "x_pts"? how many points?
  expect_equal(find_tangent_intercept(f, df, x_pts), 0)
})

## Testing different distributions
context("Test that ars works for various input functions")

# Test 2: logistic distribution
test_that("ars works for logistic distribution", {
  f <- -x + 2 * log(1 + e^x)
  df <- -1 + 2(1 - f)
  x_pts <- c(-1,1) # what value(s) go in "x_pts"? how many points?
  expect_equal(var(res_samples), var(res))
  expect_equal(mean(res_samples), mean(res))
})

# Test 3: exponential distribution
test_that("ars works for exponential distribution", {
  f <- (lambda * e) ^ -(lambda * x)
  df <- -2(lambda * e) ^ -(lambda * x)
  x_pts <- c(-1,1) # what value(s) go in "x_pts"?
  res_samples <- ars(n = 1000, f = h, dfunc = dh, T_current = c(-3, 0, 3), D_min
                     = -10, D_max = 30)
  expect_equal(var(res_samples), var(res))
  expect_equal(mean(res_samples), mean(res))
})

# Test 4: normal distribution
test_that("ars works for normal distribution", {
  f <- (-log(2 * pi * sigma ^ 2) / 2 - (y - mu) ^ 2 / (2 * sigma ^ 2))
  df <- -(y - mu)/sigma^2
  x_pts <- c(-1,1) # what value(s) go in "x_pts"?
  res_samples <- ars(n = 1000, f = h, dfunc = dh, T_current = c(-3, 0, 3), D_min
                     = -10, D_max = 30)
  expect_equal(var(res_samples), var(res))
  expect_equal(mean(res_samples), mean(res))
})


## Testing ranges
# Test 5: normal with infinite range
test_that("ars works for normal distribution with range -Inf to Inf", {
  f <- (-log(2 * pi * sigma ^ 2) / 2 - (y - mu) ^ 2 / (2 * sigma ^ 2))
  df <- -(y - mu)/sigma ^ 2
  x_pts <- c(-1,1) # what value(s) go in "x_pts"?
  res_samples <- ars(n = 1000, f = h, dfunc = dh, T_current = c(-3, 0, 3), D_min
                     = -Inf, D_max = Inf)
  expect_equal(var(res_samples), var(res))
  expect_equal(mean(res_samples), mean(res))
})

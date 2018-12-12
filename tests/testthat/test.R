# ARS Function testing
# "Test-Time" Tests

# Testing 4-6 distributions (quadratic, logistic, exponential, normal)
# Trying different ranges (-inf to inf)

## Log Concavity
context("Tests for Log Concavity")

# Test 1: quadratic function
test_that("tangent lines intercept at zero for quadratic function", {
  f <- function(x) x^2
  df <- function(x) 2*x
  x_pts <- c(-1,1)
  expect_equal(find_tangent_intercept(f, df, x_pts), c(0))
})

## Testing different distributions
context("Test that ars works for various input functions")

# Test 2: logistic distribution
test_that("ars works for logistic distribution, not specifying a domain", {
  n <- 1000
  res_samples <- ars(n, dlogis)
  true_samples <- rlogis(n)

  p_val <- ks.test(res_samples, true_samples)$p.value
  expect_gt(p_val, 0.05)

  mean_diff <- abs(mean(res_samples) - mean(true_samples))
  var_diff <- abs(var(res_samples) - var(true_samples))

  precision <- 0.3
  expect_lt(mean_diff, precision)
  expect_lt(var_diff, precision)
})

# Test 3: exponential distribution
test_that("ars works for exponential distribution, while specifying an appropriate domain", {
  n <- 1000
  res_samples <- ars(n, function(x) dexp(x, rate = 2), interval = c(0,30))
  true_samples <- rexp(n, rate = 2)

  p_val <- ks.test(res_samples, true_samples)$p.value
  expect_gt(p_val, 0.05)

  mean_diff <- abs(mean(res_samples) - mean(true_samples))
  var_diff <- abs(var(res_samples) - var(true_samples))

  precision <- 0.3
  expect_lt(mean_diff, precision)
  expect_lt(var_diff, precision)
})

# Test 4: normal distribution
test_that("ars works for normal distribution", {
  n <- 1000
  mu <- 2
  sd <- 5

  res_samples <- ars(n, function(x) dnorm(x, mean = mu, sd = sd))
  true_samples <- rnorm(n, mean = mu, sd = sd)

  p_val <- ks.test(res_samples, true_samples)$p.value
  expect_gt(p_val, 0.05)

  mean_diff <- abs(mean(res_samples) - mean(true_samples))
  relative_var_diff <- abs(var(res_samples)/sd^2 - 1)

  precision <- 0.5
  expect_lt(mean_diff, precision)
  expect_lt(relative_var_diff, precision)
})

# Test 5: normal distribution with low variance should give accurate mean
test_that("ars works for peaked normal distribution", {
  n <- 1000
  mu <- 2
  sd <- 0.01

  res_samples <- ars(n, function(x) dnorm(x, mean = mu, sd = sd), interval = c(1.8, 2.2))
  true_samples <- rnorm(n, mean = mu, sd = sd)

  p_val <- ks.test(res_samples, true_samples)$p.value
  expect_gt(p_val, 0.05)

  mean_diff <- abs(mean(res_samples) - mean(true_samples))
  var_diff <- abs(var(res_samples) - sd^2)

  precision <- 0.01
  expect_lt(mean_diff, precision)
  expect_lt(var_diff, precision)
})


# Test 6: normal distribution with unaccurate domain gives wrong result
test_that("ars works for peaked normal distribution", {
  n <- 1000

  res_samples <- ars(n, dnorm, interval = c(2, 10))
  true_samples <- rnorm(n)

  p_val <- ks.test(res_samples, true_samples)$p.value
  expect_lt(p_val, 0.05)

  mean_diff <- abs(mean(res_samples) - mean(true_samples))
  var_diff <- abs(var(res_samples) - 1)

  precision <- 0.5
  expect_gt(mean_diff, 2)
  expect_gt(var_diff, precision)
})



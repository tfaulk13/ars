# ARS Function testing
# "Test-Time" Tests


## Testing auxiliary functions
context("Tests for auxiliary functions")

# Test 1: Sanity check for find_tangent_intercept
test_that("tangent lines intercept at zero for quadratic function", {
  f <- function(x) x^2
  df <- function(x) 2*x
  x_pts <- c(-1,1)
  expect_equal(find_tangent_intercept(f, df, x_pts), c(0))
})

# Test 1.5: Check fo get_support_limit
test_that("get_support_limit does give a domain containing the mass of the distribution, for different distributions", {
  quantile <- 1E-6

  funcs <- c(dnorm, function(x) dnorm(x, mean = 10, sd = 20),
             dexp, function(x) dexp(x, rate = 3), function(x) dexp(x, rate = 100),
             dunif, function(x) dbeta(x, shape1 = 1, shape2 = 5), function(x) dgamma(x, shape = 3, rate = 2))

  for (f in funcs) {
    support <- get_support_limit(f)
    expect_lt(integrate(f, -Inf, support$D_min)$value, quantile)
    expect_lt(integrate(f, support$D_max, Inf)$value, quantile)
  }
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
test_that("ars works for peaked normal distribution with unaccurate domain gives wrong result", {
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


# Test 7: unifom distribution
test_that("ars works for uniform distribution with appropriate given starting_points", {
  n <- 1000

  res_samples <- ars(n, function(x) dunif(x), interval = c(0,1), starting_points = c(0.2, 0.8))
  true_samples <- runif(n)

  p_val <- ks.test(res_samples, true_samples)$p.value
  expect_gt(p_val, 0.05)

  mean_diff <- abs(mean(res_samples) - mean(true_samples))
  var_diff <- abs(var(res_samples) - var(true_samples))

  precision <- 0.2
  expect_lt(mean_diff, precision)
  expect_lt(var_diff, precision)
})

# Test 7: unifom distribution with wrong starting points returns error
test_that("ars returns error for uniform distribution with unappropriate given starting_points", {
  n <- 1000

  expect_error(ars(n, function(x) dunif(x), interval = c(0,1), starting_points = c(-1000, 0.8)), "Please give starting points where the density can be evaluated with a minimum precision.")
})

## Log Concavity
context("Tests for Log Concavity")

# Test 8: running ars on non log-concave functions returns error
test_that("ars returns error for non log-concave t-distribution", {
  n <- 1000

  expect_error(ars(n, function(x) dt(x, df = 2)), "Function is not log concave. Adaptive rejection sampling won't give a proper sample for this distribution.")
})


## Inputs
context("Testing inputs handling")

# Test 9: returns error for wrong inputs
test_that("ars returns error for wrong inputs", {
  n <- 1000

  # Specify a wrong domain interval returns an error
  expect_error(ars(n, dnorm, interval = c('1', '2')))
  expect_error(ars(n, dnorm, interval = c(1, 2, 3)))

  # Need to specify the number of samples required.
  expect_error(ars(dnorm))

  # Does not accept R expressions but functions as input.
  expect_error(ars(n, "x^2"))
})



## int_cir()

test_that("int_cir integrates a circular density to one", {

  f <- function(x) d_mixvmf(x = x, mu = c(0, 1), kappa = 2, p = 1)
  expect_equal(int_cir(f = f), 1, tolerance = 1e-4)

})

## int_sph()

test_that("int_sph integrates a spherical density to one", {

  f <- function(x) d_mixvmf(x = x, mu = c(0, 0, 1), kappa = 2, p = 1)
  expect_equal(int_sph(f = f), 1, tolerance = 1e-4)

})

## int_hypsph()

test_that("int_hypsph integrates a density to one", {

  f <- function(x) d_mixvmf(x = x, mu = c(0, 0, 0, 1), kappa = 2, p = 1)
  set.seed(42)
  expect_equal(int_hypsph(f = f, q = 3, M = 1e5), 1, tolerance = 5e-2)

})

test_that("int_hypsph honors na.rm (regression: was hardcoded to TRUE)", {

  set.seed(1)
  f_na <- function(x) {
    v <- d_mixvmf(x = x, mu = c(0, 0, 1), kappa = 2, p = 1)
    v[1] <- NA
    v
  }
  # With na.rm = FALSE the NA must propagate
  expect_true(is.na(int_hypsph(f = f_na, q = 2, M = 1e3, na.rm = FALSE)))
  # With na.rm = TRUE a finite value is returned
  set.seed(1)
  expect_true(is.finite(int_hypsph(f = f_na, q = 2, M = 1e3, na.rm = TRUE)))

})

test_that("Importance sampling recovers a concentrated integral", {

  # Uniform draws miss a concentrated integrand entirely, so the estimate
  # collapses to (almost) zero. Sampling half from the concentrated density and
  # half uniformly recovers it.
  q <- 10
  mu <- rbind(c(1, rep(0, q)))
  kappa <- 500
  f2 <- function(x) d_mixvmf(x = x, mu = mu, kappa = kappa, p = 1)^2

  # Exact value: int c(kappa)^2 exp(2 * kappa * x'mu) dw
  #  = c(kappa)^2 / c(2 * kappa)
  exact <- exp(2 * log_c_vmf(q = q, kappa = kappa) -
                 log_c_vmf(q = q, kappa = 2 * kappa))

  # Uniform draws understate the integral by many orders of magnitude
  set.seed(42)
  expect_lt(int_hypsph(f = f2, q = q, M = 1e4), 1e-6 * exact)

  # Importance sampling from the defensive mixture proposal recovers it
  set.seed(42)
  prop <- r_defensive_vmf(q = q, mu = mu, kappa = kappa, p = 1, M = 1e4)
  expect_equal(int_hypsph(f = f2, q = q, sim = prop$sim, f_sim = prop$f_sim),
               exact, tolerance = 5e-2)

  # The density itself must integrate to one under the same proposal
  f1 <- function(x) d_mixvmf(x = x, mu = mu, kappa = kappa, p = 1)
  expect_equal(int_hypsph(f = f1, q = q, sim = prop$sim, f_sim = prop$f_sim),
               1, tolerance = 1e-2)

})

test_that("int_hypsph importance weights are bounded by the uniform part", {

  # A spread-out integrand must still be integrated correctly under a highly
  # concentrated proposal: this is what the uniform half of the proposal buys
  q <- 10
  set.seed(42)
  prop <- r_defensive_vmf(q = q, mu = rbind(c(1, rep(0, q))), kappa = 4547,
                          p = 1, M = 1e4)
  expect_equal(int_hypsph(f = function(x) rep(1, nrow(x)), q = q,
                          sim = prop$sim, f_sim = prop$f_sim),
               w_q(q = q), tolerance = 1e-10)

})

test_that("int_hypsph accepts a sampler function and checks its dimension", {

  # sim as a function of M
  set.seed(42)
  r_unif <- function(M) normalize(matrix(rnorm(n = M * 3), nrow = M, ncol = 3))
  f <- function(x) d_mixvmf(x = x, mu = c(0, 0, 1), kappa = 2, p = 1)
  expect_equal(int_hypsph(f = f, q = 2, M = 1e4, sim = r_unif), 1,
               tolerance = 5e-2)

  # A sample on the wrong sphere is rejected
  expect_error(int_hypsph(f = f, q = 5, sim = matrix(0, nrow = 10, ncol = 3)),
               "ncol")

})

test_that("int_hypsph with a matrix sim gives common random numbers", {

  # Reusing the same draws makes the integral a deterministic function of the
  # remaining arguments, which is what optim() needs
  set.seed(42)
  sim <- normalize(matrix(rnorm(n = 1e3 * 3), nrow = 1e3, ncol = 3))
  f <- function(x, kappa) d_mixvmf(x = x, mu = c(0, 0, 1), kappa = kappa, p = 1)
  expect_identical(int_hypsph(f = f, q = 2, sim = sim, kappa = 2),
                   int_hypsph(f = f, q = 2, sim = sim, kappa = 2))

})

test_that("int_hypsph without a proposal is the plain uniform estimator", {

  # Bit for bit, adding the importance sampling arguments must perturb neither
  # the random stream nor the arithmetic of the uniform estimator
  q <- 3
  f <- function(x) d_mixvmf(x = x, mu = c(1, rep(0, q)), kappa = 2, p = 1)
  set.seed(42)
  int <- int_hypsph(f = f, q = q, M = 1e4)
  set.seed(42)
  sim <- normalize(matrix(rnorm(n = 1e4 * (q + 1)), nrow = 1e4, ncol = q + 1))
  expect_identical(int, mean(f(sim)) * w_q(q = q))

})

## int_cir(), int_sph(), int_hypsph()

test_that("integration routines also work with non-vectorized f", {

  f <- function(x) d_mixvmf(x = rbind(x), mu = c(0, 1), kappa = 1, p = 1)
  expect_equal(int_cir(f = f, f_vect = FALSE), 1, tolerance = 1e-4)

})

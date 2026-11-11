
## int_cir()

test_that("int_cir integrates a vMF density to one", {

  f <- function(x) d_mixvmf(x = x, mu = c(0, 1), kappa = 2, p = 1)
  expect_equal(int_cir(f = f), 1, tolerance = 1e-4)
  f_nonvec <- function(x) d_mixvmf(x = rbind(x), mu = c(0, 1), kappa = 1, p = 1)
  expect_equal(int_cir(f = f_nonvec, f_vect = FALSE), 1, tolerance = 1e-4)

})

## int_sph()

test_that("int_sph integrates a vMF density to one", {

  f <- function(x) d_mixvmf(x = x, mu = c(0, 0, 1), kappa = 2, p = 1)
  expect_equal(int_sph(f = f), 1, tolerance = 1e-4)

})

## int_hypsph()

test_that("int_hypsph integrates a vMF density to one", {

  f <- function(x) d_mixvmf(x = x, mu = c(0, 0, 0, 1), kappa = 2, p = 1)
  set.seed(42)
  expect_equal(int_hypsph(f = f, q = 3, M = 1e5), 1, tolerance = 5e-2)

})

test_that("int_hypsph honors na.rm", {

  f_na <- function(x) {

    v <- d_mixvmf(x = x, mu = c(0, 0, 1), kappa = 2, p = 1)
    v[1] <- NA
    v

  }
  set.seed(1)
  expect_true(is.na(int_hypsph(f = f_na, q = 2, M = 1e3, na.rm = FALSE)))
  set.seed(1)
  expect_true(is.finite(int_hypsph(f = f_na, q = 2, M = 1e3, na.rm = TRUE)))

})

test_that("int_hypsph accepts a sampler and checks the dimension of sim", {

  set.seed(42)
  r_unif <- function(M) normalize(matrix(rnorm(n = M * 3), nrow = M, ncol = 3))
  f <- function(x) d_mixvmf(x = x, mu = c(0, 0, 1), kappa = 2, p = 1)
  expect_equal(int_hypsph(f = f, q = 2, M = 1e4, sim = r_unif), 1,
               tolerance = 5e-2)
  expect_error(int_hypsph(f = f, q = 5, sim = matrix(0, nrow = 10, ncol = 3)),
               "ncol")

})

test_that("int_hypsph defaults to the plain uniform estimator", {

  q <- 3
  f <- function(x) d_mixvmf(x = x, mu = c(1, rep(0, q)), kappa = 2, p = 1)
  set.seed(42)
  int <- int_hypsph(f = f, q = q, M = 1e4)
  set.seed(42)
  sim <- normalize(matrix(rnorm(n = 1e4 * (q + 1)), nrow = 1e4, ncol = q + 1))
  expect_identical(int, mean(f(sim)) * w_q(q = q))

})

## safe_mixvmf()

test_that("safe_mixvmf recovers concentrated integrals, unlike uniform draws", {

  q <- 10
  mu <- rbind(c(1, rep(0, q)))
  kappa <- 500
  f2 <- function(x) d_mixvmf(x = x, mu = mu, kappa = kappa, p = 1)^2

  # Exact value: int c(kappa)^2 exp(2 * kappa * x'mu) dw
  #  = c(kappa)^2 / c(2 * kappa)
  exact <- exp(2 * log_c_vmf(q = q, kappa = kappa) -
                 log_c_vmf(q = q, kappa = 2 * kappa))
  set.seed(42)
  expect_lt(int_hypsph(f = f2, q = q, M = 1e4), 1e-6 * exact)
  set.seed(42)
  prop <- safe_mixvmf(q = q, mu = mu, kappa = kappa, p = 1, M = 1e4)
  expect_equal(int_hypsph(f = f2, q = q, sim = prop$sim, f_sim = prop$f_sim),
               exact, tolerance = 5e-2)
  f1 <- function(x) d_mixvmf(x = x, mu = mu, kappa = kappa, p = 1)
  expect_equal(int_hypsph(f = f1, q = q, sim = prop$sim, f_sim = prop$f_sim),
               1, tolerance = 1e-2)

})

test_that("safe_mixvmf density matches its draws, also when kappa is capped", {

  # A constant is recovered only if the proposal density matches its draws
  q <- 10
  mu <- rbind(c(1, rep(0, q)))
  f <- function(x) rep(1, nrow(x))
  set.seed(42)
  expect_silent(prop <- safe_mixvmf(q = q, mu = mu, kappa = 4547, p = 1,
                                    M = 1e4))
  expect_equal(int_hypsph(f = f, q = q, sim = prop$sim, f_sim = prop$f_sim),
               w_q(q = q), tolerance = 1e-10)
  expect_warning(prop <- safe_mixvmf(q = q, mu = mu, kappa = 500, p = 1,
                                     M = 1e4, kappa_max = 10),
                 "collapsed")
  expect_equal(int_hypsph(f = f, q = q, sim = prop$sim, f_sim = prop$f_sim),
               w_q(q = q), tolerance = 1e-2)

})

## int_dir()

test_that("int_dir dispatches to int_cir, int_sph, and int_hypsph", {

  f <- function(x) d_mixvmf(x = x, mu = c(1, rep(0, ncol(x) - 1)), kappa = 2,
                            p = 1)
  set.seed(987202226)
  sim <- normalize(matrix(rnorm(n = 4e3), nrow = 1e3, ncol = 4))
  expect_equal(int_dir(f = f, q = 1, sim = sim), int_cir(f = f))
  expect_equal(int_dir(f = f, q = 2, sim = sim), int_sph(f = f))
  expect_equal(int_dir(f = f, q = 3, sim = sim),
               int_hypsph(f = f, q = 3, sim = sim))

})

## log_sum_exp()

test_that("log_sum_exp handles edge cases", {

  expect_identical(log_sum_exp(x = numeric(0)), -Inf)
  expect_identical(log_sum_exp(x = c(-Inf, -Inf)), -Inf)
  expect_identical(log_sum_exp(x = c(1, Inf)), Inf)
  expect_equal(log_sum_exp(x = c(1000, 1000)), 1000 + log(2))
  expect_equal(log_sum_exp(x = c(-1000, -1000)), -1000 + log(2))

})

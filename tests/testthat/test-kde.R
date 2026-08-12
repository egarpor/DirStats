
## kde_dir()

test_that("kde_dir Fortran (L = NULL) and R (custom vMF kernel) paths agree", {

  set.seed(42)
  q <- 3
  samp <- rotasym::r_vMF(n = 50, mu = c(1, rep(0, q)), kappa = 2)
  x <- rbind(diag(1, q + 1), diag(-1, q + 1))
  kde_f <- kde_dir(x = x, data = samp, h = 0.5, L = NULL)
  kde_r <- kde_dir(x = x, data = samp, h = 0.5, L = function(z) exp(-z))
  expect_equal(kde_f, kde_r, tolerance = 1e-4)
  expect_length(kde_f, nrow(x))

})

test_that("kde_dir integrates to one over the sphere", {

  set.seed(1)
  samp <- rotasym::r_vMF(n = 100, mu = c(0, 0, 1), kappa = 3)
  integral <- int_sph(f = function(x) kde_dir(x = x, data = samp, h = 0.4))
  expect_equal(integral, 1, tolerance = 1e-2)

})

test_that("kde_dir accepts a single evaluation point as a vector", {

  set.seed(1)
  samp <- rotasym::r_vMF(n = 20, mu = c(0, 1), kappa = 1)
  expect_length(kde_dir(x = c(0, 1), data = samp, h = 0.5), 1)

})

test_that("kde_dir removes NAs in data", {

  set.seed(1)
  samp <- rotasym::r_vMF(n = 20, mu = c(0, 1), kappa = 1)
  samp_na <- rbind(samp, c(NA, NA))
  expect_equal(kde_dir(x = c(0, 1), data = samp_na, h = 0.5),
               kde_dir(x = c(0, 1), data = samp, h = 0.5))

})

test_that("kde_dir is finite at a bandwidth selected for concentrated data", {

  set.seed(5)
  samp <- rotasym::r_vMF(n = 50, mu = c(0, 0, 1), kappa = 1e6)
  h <- bw_dir_rot(samp)
  expect_true(is.finite(h) && h > 0)
  dens <- kde_dir(x = rbind(c(0, 0, 1), c(1, 0, 0)), data = samp, h = h)
  expect_true(all(is.finite(dens)) && all(dens >= 0))

})

## c_h()

test_that("c_h is positive and accepts vectors of bandwidths", {

  h <- c(0.2, 0.5, 1)
  expect_length(c_h(h = h, q = 2, L = NULL), 3)
  expect_true(all(c_h(h = h, q = 2, L = NULL) > 0))
  expect_equal(c_h(h = 0.5, q = 2, L = NULL),
               c_h(h = 0.5, q = 2, L = function(z) exp(-z)),
               tolerance = 1e-4)

})

test_that("c_h stays finite for small bandwidths", {

  for (q in 1:4) {
    ch <- c_h(h = c(0.005, 0.002, 0.001, 1e-4), q = q)
    expect_true(all(is.finite(ch)) && all(ch > 0))
  }
  h <- 1e-4
  expect_equal(c_h(h = h, q = 2), 1 / (2 * pi * h^2), tolerance = 1e-6)

})

## lambda_L(), b_L(), d_L()

test_that("kernel moments match between analytical and numerical vMF kernel", {

  L <- function(z) exp(-z)
  for (q in 1:3) {
    expect_equal(lambda_L(L = NULL, q = q), lambda_L(L = L, q = q),
                 tolerance = 1e-3)
    expect_equal(b_L(L = NULL, q = q), b_L(L = L, q = q), tolerance = 1e-3)
    expect_equal(d_L(L = NULL, q = q), d_L(L = L, q = q), tolerance = 1e-3)
  }

})

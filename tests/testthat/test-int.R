
## int_cir()

test_that("int_cir integrates a circular density to one", {

  f <- function(x) rotasym::d_vMF(x = x, mu = c(0, 1), kappa = 2)
  expect_equal(int_cir(f = f), 1, tolerance = 1e-4)

})

## int_sph()

test_that("int_sph integrates a spherical density to one", {

  f <- function(x) rotasym::d_vMF(x = x, mu = c(0, 0, 1), kappa = 2)
  expect_equal(int_sph(f = f), 1, tolerance = 1e-4)

})

## int_hypsph()

test_that("int_hypsph integrates a density to one", {

  f <- function(x) rotasym::d_vMF(x = x, mu = c(0, 0, 0, 1), kappa = 2)
  set.seed(42)
  expect_equal(int_hypsph(f = f, q = 3, M = 1e5), 1, tolerance = 5e-2)

})

test_that("int_hypsph honors na.rm (regression: was hardcoded to TRUE)", {

  set.seed(1)
  f_na <- function(x) {
    v <- rotasym::d_vMF(x = x, mu = c(0, 0, 1), kappa = 2)
    v[1] <- NA
    v
  }
  # With na.rm = FALSE the NA must propagate
  expect_true(is.na(int_hypsph(f = f_na, q = 2, M = 1e3, na.rm = FALSE)))
  # With na.rm = TRUE a finite value is returned
  set.seed(1)
  expect_true(is.finite(int_hypsph(f = f_na, q = 2, M = 1e3, na.rm = TRUE)))

})

## int_cir(), int_sph(), int_hypsph()

test_that("integration routines also work with non-vectorized f", {

  f <- function(x) rotasym::d_vMF(x = rbind(x), mu = c(0, 1), kappa = 1)
  expect_equal(int_cir(f = f, f_vect = FALSE), 1, tolerance = 1e-4)

})

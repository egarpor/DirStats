
## log_besselI_scaled()

test_that("log_besselI_scaled matches base besselI in the standard regime", {

  x <- exp(seq(log(1e-3), log(5e3), length.out = 200))
  for (nu in c(0, 0.5, 1.5, 10, 24.5, 40)) {
    ref <- log(besselI(x, nu, expon.scaled = TRUE))
    expect_equal(log_besselI_scaled(nu = nu, x = x, spline = FALSE), ref,
                 tolerance = 1e-10)
  }

})

test_that("log_besselI_scaled spline path matches the direct evaluation", {

  x <- exp(seq(log(1e-4), log(1e4), length.out = 500))
  for (nu in c(0, 0.5, 1.5, 24.5, 40, 50)) {
    spl <- log_besselI_scaled(nu = nu, x = x, spline = TRUE)
    dir <- log_besselI_scaled(nu = nu, x = x, spline = FALSE)
    expect_equal(spl, dir, tolerance = 1e-6)
  }

})

test_that("log_besselI_scaled joins smoothly to the large-order asymptotic", {

  # At the table cutoff (nu = 50) the spline/table value and the large-order
  # asymptotic that takes over for nu > 50 must agree closely
  x <- exp(seq(log(1e-3), log(1e4), length.out = 300))
  spl <- log_besselI_scaled(nu = 50, x = x, spline = TRUE)
  asy <- Bessel::besselI.nuAsym(x = x, nu = 50, k.max = 5,
                                expon.scaled = TRUE, log = TRUE)
  expect_lt(max(abs(spl - asy)), 1e-8)

})

test_that("log_besselI_scaled stays finite for large order and argument", {

  # Where base besselI() overflows/underflows or warns, the asymptotics keep
  # the result finite and warning-free
  expect_silent(v1 <- log_besselI_scaled(nu = 800, x = 20))
  expect_true(is.finite(v1))
  expect_silent(v2 <- log_besselI_scaled(nu = 1.5, x = 1e6))
  expect_true(is.finite(v2))
  # Very large order and argument together
  expect_true(is.finite(log_besselI_scaled(nu = 500, x = 1e5)))

})

test_that("log_besselI_scaled vectorizes and propagates NAs", {

  expect_length(log_besselI_scaled(nu = c(0.5, 1.5, 2.5), x = c(1, 2, 3)), 3)
  expect_length(log_besselI_scaled(nu = 0.5, x = 1:4), 4)
  out <- log_besselI_scaled(nu = 0.5, x = c(1, NA, 3))
  expect_true(is.na(out[2]) && all(is.finite(out[-2])))
  expect_error(log_besselI_scaled(nu = c(1, 2), x = c(1, 2, 3)), "scalar")
  # Zero-length input gives zero-length output, as besselI() does
  expect_length(log_besselI_scaled(nu = 0.5, x = numeric(0)), 0)

})

test_that("The fast path of log_besselI_scaled() matches the general path", {

  # An NA in x forces the general path, which must give identical values for
  # the non-NA entries. besselI() warns where it underflows (nu = 50 and
  # x = 1e-5), which both paths fill with the small-argument series
  x <- c(1e-5, 1e-4, 0.5, 3, 40, 900, 9999)
  for (spline in c(TRUE, FALSE)) {

    for (nu in list(0, 1.5, 7, 50, -0.5, 0:6, rep(1.5, 7))) {

      fast <- suppressWarnings(log_besselI_scaled(nu = nu, x = x,
                                                  spline = spline))
      nu_na <- if (length(nu) > 1) c(nu, nu[1]) else nu
      general <- suppressWarnings(log_besselI_scaled(nu = nu_na, x = c(x, NA),
                                                     spline = spline))
      expect_identical(fast, general[seq_along(x)])

    }

  }

  # Matrix input gives a plain vector in both paths
  expect_identical(log_besselI_scaled(nu = 1, x = matrix(1:4, 2),
                                      spline = FALSE),
                   log_besselI_scaled(nu = 1, x = 1:4, spline = FALSE))

})

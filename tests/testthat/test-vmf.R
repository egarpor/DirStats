test_that("mu_ml recovers the mean direction", {

  set.seed(42)
  mu <- normalize(c(1, 2, -1))
  samp <- rotasym::r_vMF(n = 1000, mu = mu, kappa = 20)
  expect_equal(mu_ml(samp), mu, tolerance = 5e-2)
  expect_equal(norm2(mu_ml(samp)), 1)

})


test_that("kappa_ml recovers the concentration", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 2000, mu = c(0, 0, 1), kappa = 5)
  expect_equal(kappa_ml(samp), 5, tolerance = 0.5)
  # A search interval that does not bracket the root cannot be solved: the
  # min_kappa / max_kappa arguments are passed to uniroot()
  expect_error(kappa_ml(samp, max_kappa = 0.5), "opposite sign")

})


test_that("d_mixvmf integrates to one on the circle", {

  x <- to_cir(seq(0, 2 * pi, l = 500))
  mu <- rbind(c(-1, 0), c(0, 1), c(1, 0))
  dens <- function(y) d_mixvmf(x = y, mu = mu, kappa = 1:3, p = c(0.5, 0.2, 0.3))
  expect_equal(int_cir(f = dens), 1, tolerance = 1e-3)

})


test_that("d_mixvmf with a single component equals rotasym::d_vMF", {

  x <- to_cir(seq(0, 2 * pi, l = 7))
  # d_mixvmf returns a plain vector; rotasym::d_vMF returns a column matrix
  expect_equal(d_mixvmf(x = x, mu = rbind(c(1, 0)), kappa = 2, p = 1),
               rotasym::d_vMF(x = x, mu = c(1, 0), kappa = 2),
               ignore_attr = TRUE)

})


test_that("d_mixvmf checks the mixture argument dimensions", {

  x <- to_cir(c(0, 1))
  expect_error(d_mixvmf(x = x, mu = rbind(c(1, 0), c(0, 1)),
                        kappa = 1, p = c(0.5, 0.5)), "mixture")

})


test_that("d_mixvmf handles a single evaluation point, norm and vector inputs", {

  mu <- rbind(c(-1, 0), c(0, 1))
  # x given as a plain vector (single point), with norm = TRUE, and a
  # multi-component mixture (so the density reshape branch is exercised)
  d_vec <- d_mixvmf(x = c(0, 2), mu = mu, kappa = c(1, 2),
                    p = c(0.5, 0.5), norm = TRUE)
  expect_length(d_vec, 1)
  # Must match the matrix-input evaluation at the normalized point
  d_mat <- d_mixvmf(x = normalize(rbind(c(0, 2))), mu = mu, kappa = c(1, 2),
                    p = c(0.5, 0.5))
  expect_equal(d_vec, d_mat)

})


test_that("d_mixvmf accepts a single component given as a vector mu", {

  x <- to_cir(c(0, 1, 2))
  expect_equal(d_mixvmf(x = x, mu = c(1, 0), kappa = 2, p = 1),
               rotasym::d_vMF(x = x, mu = c(1, 0), kappa = 2),
               ignore_attr = TRUE)

})

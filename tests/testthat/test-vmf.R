
## mu_ml()

test_that("mu_ml recovers the mean direction", {

  set.seed(42)
  mu <- normalize(c(1, 2, -1))
  samp <- r_proj_norm(n = 1000, mu = mu, kappa = 20)
  expect_equal(mu_ml(samp), mu, tolerance = 5e-2)
  expect_equal(norm2(mu_ml(samp)), 1)

})

## kappa_ml()

test_that("kappa_ml recovers the concentration for large kappa and dimension", {

  for (q in c(1:2, 100)) {
    for (true_kappa in c(2, 50, 300)) {

      set.seed(42)
      samp <- movMF::rmovMF(n = 1e4, theta = rbind(c(true_kappa, rep(0, q))))
      k_hat <- kappa_ml(samp)
      expect_true(is.finite(k_hat) && k_hat > 0)
      expect_equal(k_hat, true_kappa, tolerance = 1e-1)

    }
  }

})

test_that("kappa_ml extends the search beyond a max_kappa below kappa", {

  set.seed(42)
  samp <- movMF::rmovMF(n = 2000, theta = rbind(c(0, 0, 5)))
  expect_equal(kappa_ml(samp), 5, tolerance = 0.5)
  expect_equal(kappa_ml(samp, max_kappa = 0.5), kappa_ml(samp),
               tolerance = 1e-3)

})

test_that("kappa_ml accepts a pre-computed resultant length", {

  set.seed(42)
  samp <- movMF::rmovMF(n = 2000, theta = rbind(c(0, 0, 5)))
  expect_equal(kappa_ml(samp, R = norm2(colMeans(samp))), kappa_ml(samp))

})

test_that("kappa_ml handles degenerate resultant lengths", {

  ident <- matrix(c(0, 0, 1), nrow = 5, ncol = 3, byrow = TRUE)
  antipodal <- rbind(c(0, 0, 1), c(0, 0, -1))
  expect_equal(kappa_ml(ident), Inf)
  expect_equal(kappa_ml(antipodal), 0)

})

## d_mixvmf()

test_that("d_mixvmf integrates to one on the circle", {

  x <- to_cir(seq(0, 2 * pi, l = 500))
  mu <- rbind(c(-1, 0), c(0, 1), c(1, 0))
  dens <- function(y) d_mixvmf(x = y, mu = mu, kappa = 1:3,
                               p = c(0.5, 0.2, 0.3))
  expect_equal(int_cir(f = dens), 1, tolerance = 1e-3)

})

test_that("d_mixvmf matches the von Mises density for matrix or vector mu", {

  x <- to_cir(seq(0, 2 * pi, l = 7))
  d_vm <- exp(2 * x[, 1]) / (2 * pi * besselI(x = 2, nu = 0))
  expect_equal(d_mixvmf(x = x, mu = rbind(c(1, 0)), kappa = 2, p = 1), d_vm,
               ignore_attr = TRUE)
  expect_equal(d_mixvmf(x = x, mu = c(1, 0), kappa = 2, p = 1), d_vm,
               ignore_attr = TRUE)

})

test_that("d_mixvmf normalizes a vector x if norm = TRUE", {

  mu <- rbind(c(-1, 0), c(0, 1))
  d_vec <- d_mixvmf(x = c(0, 2), mu = mu, kappa = c(1, 2),
                    p = c(0.5, 0.5), norm = TRUE)
  expect_length(d_vec, 1)
  d_mat <- d_mixvmf(x = normalize(rbind(c(0, 2))), mu = mu, kappa = c(1, 2),
                    p = c(0.5, 0.5))
  expect_equal(d_vec, d_mat)

})

test_that("d_mixvmf stays finite for large concentrations", {

  for (q in 1:3) {

    mu <- rbind(c(rep(0, q), 1))
    x <- rbind(c(rep(0, q), 1), c(1, rep(0, q)))
    dens <- sapply(c(1e4, 1e5, 1e6, 1e8), function(k)
      d_mixvmf(x = x, mu = mu, kappa = k, p = 1))
    expect_true(all(is.finite(dens)))

    # The density at the mean direction grows with the concentration, and
    # vanishes away from it.
    expect_true(all(diff(dens[1, ]) > 0))
    expect_equal(dens[2, ], rep(0, 4))

  }

})

test_that("d_mixvmf checks the mixture argument dimensions", {

  x <- to_cir(c(0, 1))
  expect_error(d_mixvmf(x = x, mu = rbind(c(1, 0), c(0, 1)),
                        kappa = 1, p = c(0.5, 0.5)), "mixture")

})

## log_c_vmf()

test_that("log_c_vmf normalizes the density, also vectorized in q", {

  # 1 = c_q(kappa) * w_q *
  #     int_0^pi exp(kappa * cos(th)) * sin(th)^(q - 1) dth
  for (q in c(1, 2, 4, 9)) {
    for (kappa in c(0, 0.5, 5, 50, 300)) {

      int <- integrate(f = function(th) {
        exp(kappa * (cos(th) - 1)) * sin(th)^(q - 1)
      }, lower = 0, upper = pi, rel.tol = 1e-10)$value
      expect_equal(log_c_vmf(q = q, kappa = kappa) +
                     w_q(q = q - 1, log = TRUE) + kappa + log(int),
                   0, tolerance = 1e-8)

    }
  }

  # A vector q, with some zero concentrations
  q <- 1:3
  kappa <- c(0, 1, 0)
  expect_equal(log_c_vmf(q = q, kappa = kappa),
               sapply(seq_along(q), function(i)
                 log_c_vmf(q = q[i], kappa = kappa[i])))
  expect_equal(log_c_vmf(q = q, kappa = 0), -w_q(q = q, log = TRUE))

})

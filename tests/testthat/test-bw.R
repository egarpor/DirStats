test_that("bw_dir_rot returns a positive bandwidth for q = 1, 2, 3", {

  set.seed(42)
  for (q in 1:3) {
    samp <- rotasym::r_vMF(n = 50, mu = c(1, rep(0, q)), kappa = 2)
    h <- bw_dir_rot(samp)
    expect_length(h, 1)
    expect_true(is.finite(h) && h > 0)
  }

})


test_that("bw_dir_lcv Fortran and R (custom vMF kernel) paths agree on a grid", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 30, mu = c(1, 0, 0), kappa = 2)
  h_grid <- exp(seq(log(0.1), log(1.2), l = 15))

  cv_f <- bw_dir_lcv(data = samp, h_grid = h_grid, optim = FALSE)
  cv_r <- bw_dir_lcv(data = samp, h_grid = h_grid, optim = FALSE,
                     L = function(z) exp(-z))

  expect_equal(cv_f$CV_grid, cv_r$CV_grid, tolerance = 1e-3)
  expect_equal(cv_f$h_opt, cv_r$h_opt)

})


test_that("bw_dir_lcv optim and grid search give consistent bandwidths", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 40, mu = c(1, 0, 0), kappa = 2)
  h_opt_optim <- bw_dir_lcv(data = samp, optim = TRUE)$h_opt
  h_opt_grid <- bw_dir_lcv(data = samp, optim = FALSE)$h_opt
  expect_equal(h_opt_optim, h_opt_grid, tolerance = 0.1)

})


test_that("bw_dir_lscv Fortran and R paths agree on a grid", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 25, mu = c(1, 0, 0), kappa = 2)
  h_grid <- exp(seq(log(0.2), log(1.2), l = 12))

  lscv_f <- bw_dir_lscv(data = samp, h_grid = h_grid, optim = FALSE)
  lscv_r <- bw_dir_lscv(data = samp, h_grid = h_grid, optim = FALSE,
                        R_code = TRUE)

  expect_equal(lscv_f$CV_grid, lscv_r$CV_grid, tolerance = 1e-4)

})


test_that("bw selectors warn when the optimum is at the grid boundary", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 20, mu = c(1, 0, 0), kappa = 2)
  # A tiny grid of large bandwidths forces the optimum to a boundary
  expect_warning(bw_dir_lcv(data = samp, h_grid = c(1.4, 1.5), optim = FALSE),
                 "extreme")

})


test_that("bic_vmf_mix returns the documented structure", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 60, mu = c(1, 0, 0), kappa = 3)
  fit <- bic_vmf_mix(data = samp)
  expect_named(fit, c("best_fit", "fits", "BICs"))
  expect_named(fit$best_fit, c("mu_hat", "kappa_hat", "p_hat"))
  expect_equal(sum(fit$best_fit$p_hat), 1, tolerance = 1e-6)

})


test_that("bic_vmf_mix supports BIC, AIC and AICc criteria (AICc regression)", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 60, mu = c(1, 0, 0), kappa = 3)
  for (crit in c("BIC", "AIC", "AICc")) {
    fit <- bic_vmf_mix(data = samp, crit = crit, iterative = FALSE)
    expect_true(is.list(fit))
    expect_true(is.numeric(fit$best_fit$kappa_hat))
  }

})


test_that("R_Psi_mixvmf agrees between exact and numeric routes", {

  # Single component uses the exact expression; wrapping the same component
  # twice with equal weights must reproduce it (numeric integration route)
  mu <- rbind(c(0, 0, 1))
  R_exact <- R_Psi_mixvmf(q = 2, mu = mu, kappa = 4, p = 1)
  R_numer <- R_Psi_mixvmf(q = 2, mu = rbind(c(0, 0, 1), c(0, 0, 1)),
                          kappa = c(4, 4), p = c(0.5, 0.5))
  expect_equal(R_exact, R_numer, tolerance = 1e-2)

})


test_that("R_Psi_mixvmf checks argument dimensions", {

  expect_error(R_Psi_mixvmf(q = 2, mu = rbind(c(0, 0, 1)), kappa = c(1, 2),
                            p = c(0.5, 0.5)), "dimension")

})


test_that("bw_dir_ami returns a positive bandwidth", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 50, mu = c(1, 0, 0), kappa = 3)
  fit_mix <- bic_vmf_mix(data = samp)
  h <- bw_dir_ami(samp, fit_mix = fit_mix)
  expect_length(h, 1)
  expect_true(is.finite(h) && h > 0)

})


test_that("bw_dir_emi returns a positive bandwidth", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 40, mu = c(1, 0, 0), kappa = 3)
  fit_mix <- bic_vmf_mix(data = samp)
  emi <- bw_dir_emi(samp, fit_mix = fit_mix, optim = TRUE)
  expect_true(is.finite(emi$h_opt) && emi$h_opt > 0)

})

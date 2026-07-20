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


# ---------------------------------------------------------------------------
# Additional coverage for the alternative branches of the CV selectors
# ---------------------------------------------------------------------------

test_that("bw_dir_lscv optimizes over the Fortran vMF path", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 25, mu = c(1, 0, 0), kappa = 2)
  lscv <- bw_dir_lscv(data = samp, optim = TRUE)
  expect_null(lscv$h_grid)
  expect_true(is.finite(lscv$h_opt) && lscv$h_opt > 0)

})


test_that("bw_dir_lscv with a custom kernel works on the circle (q = 1)", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 20, mu = c(1, 0), kappa = 2)
  h_grid <- exp(seq(log(0.2), log(1.2), l = 8))
  # A boundary warning may fire depending on the sample; not what is tested here
  lscv <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                       optim = FALSE, L = function(z) exp(-z)))
  expect_length(lscv$CV_grid, length(h_grid))
  expect_true(all(is.finite(lscv$CV_grid)))

})


test_that("bw_dir_lscv with a custom kernel works on the sphere (q = 2)", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 20, mu = c(1, 0, 0), kappa = 2)
  h_grid <- exp(seq(log(0.2), log(1.2), l = 8))
  lscv <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                       optim = FALSE, L = function(z) exp(-z)))
  expect_length(lscv$CV_grid, length(h_grid))
  expect_true(all(is.finite(lscv$CV_grid)))

})


test_that("bw_dir_lscv with a custom kernel works on S^3 (q = 3)", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 20, mu = c(1, 0, 0, 0), kappa = 2)
  h_grid <- exp(seq(log(0.3), log(1.2), l = 6))
  lscv <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                       optim = FALSE, L = function(z) exp(-z)))
  expect_length(lscv$CV_grid, length(h_grid))
  expect_true(all(is.finite(lscv$CV_grid)))

})


test_that("the CV selectors honor plot_it (grid search)", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 20, mu = c(1, 0, 0), kappa = 2)
  h_grid <- exp(seq(log(0.1), log(1.2), l = 10))
  # Plots are routed to a null device (see setup.R); a boundary warning may or
  # may not fire depending on the sample, so it is suppressed here
  lcv <- suppressWarnings(bw_dir_lcv(data = samp, h_grid = h_grid,
                                     optim = FALSE, plot_it = TRUE))
  lscv <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                       optim = FALSE, plot_it = TRUE,
                                       L = function(z) exp(-z)))
  expect_true(is.finite(lcv$h_opt))
  expect_true(is.finite(lscv$h_opt))

})


test_that("the CV selectors warn at the boundary of a short grid", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 20, mu = c(1, 0, 0), kappa = 2)
  expect_warning(bw_dir_lscv(data = samp, h_grid = c(1.4, 1.5), optim = FALSE),
                 "extreme")

})


# ---------------------------------------------------------------------------
# Additional coverage for bic_vmf_mix branches
# ---------------------------------------------------------------------------

test_that("bic_vmf_mix runs with verbose output", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 60, mu = c(1, 0, 0), kappa = 3)
  out <- capture.output(fit <- bic_vmf_mix(data = samp, verbose = TRUE))
  expect_true(any(grepl("Fitting mixtures", out)))
  expect_named(fit, c("best_fit", "fits", "BICs"))

})


test_that("bic_vmf_mix runs the iterative extension of the search", {

  # A well-separated three-component sample, with a low M_bound and M_neig = 3
  # so that (M_bound - M_neig + 1) <= M_opt and the repeat{} loop is entered
  set.seed(42)
  q <- 2
  samp <- rbind(
    rotasym::r_vMF(n = 30, mu = c(rep(0, q), 1), kappa = 8),
    rotasym::r_vMF(n = 30, mu = c(rep(0, q), -1), kappa = 8),
    rotasym::r_vMF(n = 30, mu = c(1, rep(0, q)), kappa = 8))
  fit <- bic_vmf_mix(data = samp, M_bound = 2, M_neig = 3, iterative = TRUE)
  expect_named(fit, c("best_fit", "fits", "BICs"))
  expect_true(length(fit$BICs) > 2)
  expect_gte(length(fit$best_fit$kappa_hat), 1)

})


test_that("bic_vmf_mix honors plot_it", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 60, mu = c(1, 0, 0), kappa = 3)
  # Plot goes to the null device (see setup.R)
  fit <- bic_vmf_mix(data = samp, plot_it = TRUE)
  expect_named(fit, c("best_fit", "fits", "BICs"))

})


test_that("bic_vmf_mix falls back to a single vMF when all fits are rejected", {

  # A tiny kappa_max rejects every fitted mixture, so which.min(BICs) is empty
  # and the maximum-likelihood single-vMF fallback is used (iterative = FALSE
  # avoids comparing against the empty M_opt)
  set.seed(42)
  samp <- rotasym::r_vMF(n = 40, mu = c(1, 0, 0), kappa = 3)
  fit <- bic_vmf_mix(data = samp, kappa_max = 1e-3, iterative = FALSE)
  expect_equal(fit$best_fit$p_hat, 1)
  expect_equal(norm2(fit$best_fit$mu_hat), 1)
  expect_true(is.finite(fit$best_fit$kappa_hat))

})


# ---------------------------------------------------------------------------
# Additional coverage for R_Psi_mixvmf, bw_dir_ami and bw_dir_emi branches
# ---------------------------------------------------------------------------

test_that("R_Psi_mixvmf numeric route matches the exact one on the circle", {

  # A single component uses the exact expression; duplicating it with equal
  # weights must reproduce it through the numeric (int_cir) route for q = 1
  R_exact <- R_Psi_mixvmf(q = 1, mu = rbind(c(1, 0)), kappa = 3, p = 1)
  R_numer <- R_Psi_mixvmf(q = 1, mu = rbind(c(1, 0), c(1, 0)),
                          kappa = c(3, 3), p = c(0.5, 0.5))
  expect_equal(R_exact, R_numer, tolerance = 1e-2)

})


test_that("R_Psi_mixvmf numeric route works on S^3 (q = 3)", {

  set.seed(42)
  R_exact <- R_Psi_mixvmf(q = 3, mu = rbind(c(1, 0, 0, 0)), kappa = 3, p = 1)
  R_numer <- R_Psi_mixvmf(q = 3, mu = rbind(c(1, 0, 0, 0), c(1, 0, 0, 0)),
                          kappa = c(3, 3), p = c(0.5, 0.5))
  expect_equal(R_exact, R_numer, tolerance = 5e-2)

})


test_that("bw_dir_ami fits the mixture internally when fit_mix is NULL", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 50, mu = c(1, 0, 0), kappa = 3)
  h <- bw_dir_ami(samp)
  expect_true(is.finite(h) && h > 0)

})


test_that("bw_dir_emi fits the mixture internally and integrates on S^1", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 40, mu = c(1, 0), kappa = 3)
  emi <- bw_dir_emi(samp, optim = TRUE)
  expect_true(is.finite(emi$h_opt) && emi$h_opt > 0)

})


test_that("bw_dir_emi grid search works with a single-component mixture", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 40, mu = c(0, 0, 1), kappa = 3)
  fit_mix <- list(best_fit = list(mu_hat = rbind(c(0, 0, 1)),
                                  kappa_hat = 3, p_hat = 1))
  h_grid <- exp(seq(log(0.1), log(1.2), l = 8))
  emi <- suppressWarnings(bw_dir_emi(samp, fit_mix = fit_mix, optim = FALSE,
                                     plot_it = TRUE, h_grid = h_grid))
  expect_length(emi$MISE_grid, length(h_grid))
  expect_true(is.finite(emi$h_opt))

})


test_that("bw_dir_emi integrates on S^3 (q = 3)", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 30, mu = c(1, 0, 0, 0), kappa = 3)
  fit_mix <- list(best_fit = list(mu_hat = rbind(c(1, 0, 0, 0)),
                                  kappa_hat = 3, p_hat = 1))
  emi <- bw_dir_emi(samp, fit_mix = fit_mix, optim = TRUE)
  expect_true(is.finite(emi$h_opt) && emi$h_opt > 0)

})


test_that("bw_dir_emi warns at the boundary of a short grid", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 30, mu = c(0, 0, 1), kappa = 3)
  fit_mix <- list(best_fit = list(mu_hat = rbind(c(0, 0, 1)),
                                  kappa_hat = 3, p_hat = 1))
  expect_warning(bw_dir_emi(samp, fit_mix = fit_mix, optim = FALSE,
                            plot_it = FALSE, h_grid = c(1.4, 1.5)),
                 "extreme")

})

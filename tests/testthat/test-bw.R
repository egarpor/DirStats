
## bw_dir_rot()

test_that("bw_dir_rot returns a positive bandwidth for q = 1, 2, 3", {

  set.seed(42)
  for (q in 1:3) {
    samp <- rotasym::r_vMF(n = 50, mu = c(1, rep(0, q)), kappa = 2)
    h <- bw_dir_rot(samp)
    expect_length(h, 1)
    expect_true(is.finite(h) && h > 0)
  }

})

test_that("bw_dir_rot q = 1 is unchanged by the scaled-Bessel rewrite", {

  # Reference values produced by the previous implementation
  set.seed(123)
  samp <- rotasym::r_vMF(n = 100, mu = c(1, 0), kappa = 5)
  h <- bw_dir_rot(samp)
  expect_equal(h, 0.19110687814973268, tolerance = 1e-6)
  expect_true(is.finite(h) && h > 0)

})

test_that("bw_dir_rot log-space formula reproduces the closed-form values", {

  # Reference values produced by the previous implementation
  ref <- c(0.22344720788485600, 0.24192194154685706,
           0.23699316703400239, 0.26285215925296990)
  for (q in 1:4) {
    set.seed(200 + q)
    samp <- rotasym::r_vMF(n = 80, mu = c(1, rep(0, q)), kappa = 4)
    expect_equal(bw_dir_rot(samp), ref[q], tolerance = 1e-8)
  }

})

test_that("bw_dir_rot stays finite for very large concentrations", {

  set.seed(5)
  for (true_kappa in c(1e5, 1e6)) {
    for (q in 1:3) {
      samp <- rotasym::r_vMF(n = 50, mu = c(1, rep(0, q)), kappa = true_kappa)
      h <- bw_dir_rot(samp)
      expect_true(is.finite(h) && h > 0)
    }
  }

})

test_that("bw_dir_rot warns and returns NA for degenerate concentrations", {

  ident <- matrix(c(0, 0, 1), nrow = 5, ncol = 3, byrow = TRUE)
  expect_warning(h <- bw_dir_rot(ident), "[Dd]egenerate")
  expect_true(is.na(h))

})

## bw_dir_lcv()

# Naive reference for the von Mises-Fisher LCV
lcv_naive <- function(data, h) {

  n <- nrow(data)
  p <- ncol(data)
  f <- sapply(h, function(hh) sapply(seq_len(n), function(i) {
    rotasym::c_vMF(p = p, kappa = 1 / hh^2) / (n - 1) *
      sum(exp(data[-i, ] %*% data[i, ] / hh^2))
  }))
  return(-colMeans(log(f)))

}

test_that("bw_dir_lcv and the naive likelihood CV loss agree", {

  h_grid <- exp(seq(log(0.4), log(1.2), l = 5))
  for (q in 1:3) {

    set.seed(7)
    samp <- rotasym::r_vMF(n = 12, mu = c(1, rep(0, q)), kappa = 2)
    lcv <- suppressWarnings(bw_dir_lcv(data = samp, h_grid = h_grid,
                                       optim = FALSE))
    expect_equal(lcv$CV_grid, lcv_naive(data = samp, h = h_grid),
                 tolerance = 1e-8)

  }

})

test_that("bw_dir_lcv von Mises-Fisher and custom kernel paths agree", {

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

test_that("bw_dir_lcv warns when the optimum is at the grid boundary", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 20, mu = c(1, 0, 0), kappa = 2)
    expect_warning(bw_dir_lcv(data = samp, h_grid = c(1.4, 1.5), optim = FALSE),
                 "extreme")

})


## bw_dir_lscv()

# Naive reference for the von Mises-Fisher LSCV
lscv_naive <- function(data, h) {

  n <- nrow(data)
  p <- ncol(data)
  return(sapply(h, function(hh) {

    ck <- rotasym::c_vMF(p = p, kappa = 1 / hh^2)
    ck2 <- rotasym::c_vMF(p = p, kappa = 2 / hh^2)
    s1 <- s2 <- 0
    for (i in 2:n) {
      for (j in 1:(i - 1)) {

        s1 <- s1 + exp(sum(data[i, ] * data[j, ]) / hh^2)
        s2 <- s2 + 1 / rotasym::c_vMF(p = p,
                                      kappa = norm2(data[i, ] + data[j, ]) /
                                        hh^2)

      }
    }
    -(4 * ck / (n * (n - 1)) * s1 - ck^2 / (n * ck2) - 2 * ck^2 / n^2 * s2)

  }))

}

test_that("bw_dir_lscv and the naive least squares CV loss agree", {

  h_grid <- exp(seq(log(0.4), log(1.2), l = 5))
  for (q in 1:3) {

    set.seed(7)
    samp <- rotasym::r_vMF(n = 12, mu = c(1, rep(0, q)), kappa = 2)
    lscv <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                         optim = FALSE))
    expect_equal(lscv$CV_grid, lscv_naive(data = samp, h = h_grid),
                 tolerance = 1e-7)

  }

})

test_that("bw_dir_lscv optimizes over the von Mises-Fisher path", {

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

test_that("the CV selectors warn at the boundary of a short grid", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 20, mu = c(1, 0, 0), kappa = 2)
  expect_warning(bw_dir_lscv(data = samp, h_grid = c(1.4, 1.5), optim = FALSE),
                 "extreme")

})

test_that("bw_dir_lscv Fortran and R paths agree below the RIBESL crossover", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 25, mu = c(1, 0, 0), kappa = 2)
  h_grid <- exp(seq(log(0.005), log(0.5), l = 15))

  lscv_f <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                         optim = FALSE))
  lscv_r <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                         optim = FALSE, R_code = TRUE))
  expect_equal(lscv_f$CV_grid, lscv_r$CV_grid, tolerance = 1e-6)

})

## bic_vmf_mix()

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

test_that("bic_vmf_mix runs with verbose output", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 60, mu = c(1, 0, 0), kappa = 3)
  out <- capture.output(fit <- bic_vmf_mix(data = samp, verbose = TRUE))
  expect_true(any(grepl("Fitting mixtures", out)))
  expect_named(fit, c("best_fit", "fits", "BICs"))

})

test_that("bic_vmf_mix runs the iterative extension of the search", {

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
  fit <- bic_vmf_mix(data = samp, plot_it = TRUE)
  expect_named(fit, c("best_fit", "fits", "BICs"))

})

## R_Psi_mixvmf()

test_that("R_Psi_mixvmf agrees between exact and numeric routes", {

  # Single component uses the exact expression; wrapping the same component
  # twice with equal weights must reproduce it (numeric integration route)
  mu <- rbind(c(0, 0, 1))
  R_exact <- R_Psi_mixvmf(q = 2, mu = mu, kappa = 4, p = 1)
  R_numer <- R_Psi_mixvmf(q = 2, mu = rbind(c(0, 0, 1), c(0, 0, 1)),
                          kappa = c(4, 4), p = c(0.5, 0.5))
  expect_equal(R_exact, R_numer, tolerance = 1e-2)

})

test_that("R_Psi_mixvmf stays finite for large concentrations", {

  for (q in 1:3) {
    mu <- rbind(c(rep(0, q), 1))
    R <- sapply(c(1e4, 1e5, 1e6, 1e8), function(k)
      R_Psi_mixvmf(q = q, mu = mu, kappa = k, p = 1))
    expect_true(all(is.finite(R)) && all(R > 0))
    # The curvature must be increasing in the concentration
    expect_true(all(diff(R) > 0))
  }

})

test_that("R_Psi_mixvmf checks argument dimensions", {

  expect_error(R_Psi_mixvmf(q = 2, mu = rbind(c(0, 0, 1)), kappa = c(1, 2),
                            p = c(0.5, 0.5)), "dimension")

})

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

## bw_dir_ami()

test_that("bw_dir_ami returns a positive bandwidth", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 50, mu = c(1, 0, 0), kappa = 3)
  fit_mix <- bic_vmf_mix(data = samp)
  h <- bw_dir_ami(samp, fit_mix = fit_mix)
  expect_length(h, 1)
  expect_true(is.finite(h) && h > 0)

})

test_that("bw_dir_ami fits the mixture internally when fit_mix is NULL", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 50, mu = c(1, 0, 0), kappa = 3)
  h <- bw_dir_ami(samp)
  expect_true(is.finite(h) && h > 0)

})

## bw_dir_emi()

test_that("bw_dir_emi returns a positive bandwidth", {

  set.seed(42)
  samp <- rotasym::r_vMF(n = 40, mu = c(1, 0, 0), kappa = 3)
  fit_mix <- bic_vmf_mix(data = samp)
  emi <- bw_dir_emi(samp, fit_mix = fit_mix, optim = TRUE)
  expect_true(is.finite(emi$h_opt) && emi$h_opt > 0)

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

test_that("bw_dir_emi MISE is finite on a grid extending below h = 0.05", {

  h_grid <- exp(seq(log(0.03), log(1.5), l = 15))
  for (q in 1:3) {
    set.seed(7)
    samp <- rotasym::r_vMF(n = 30, mu = c(1, rep(0, q)), kappa = 3)
    fit_mix <- list(best_fit = list(mu_hat = rbind(mu_ml(samp)),
                                    kappa_hat = 3, p_hat = 1))
    set.seed(7)
    mise <- suppressWarnings(
      bw_dir_emi(samp, fit_mix = fit_mix, optim = FALSE, plot_it = FALSE,
                 h_grid = h_grid)$MISE_grid)
    expect_length(mise, length(h_grid))
    expect_true(all(is.finite(mise)))
    expect_true(all(mise > 0))
  }

})

## CV selectors (bw_dir_rot(), bw_dir_lcv(), bw_dir_lscv())

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

test_that("The CV losses stay finite for arbitrarily small bandwidths", {

  h_grid <- c(sqrt(.Machine$double.eps), 1e-6, 1e-3, 0.01, 0.05, 0.5)
  for (q in 1:3) {

    set.seed(7)
    samp <- rotasym::r_vMF(n = 20, mu = c(1, rep(0, q)), kappa = 3)
    lcv <- suppressWarnings(bw_dir_lcv(data = samp, h_grid = h_grid,
                                       optim = FALSE)$CV_grid)
    lscv <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                         optim = FALSE)$CV_grid)
    expect_true(all(is.finite(lcv)))
    expect_true(all(is.finite(lscv)))

  }

})

test_that("The CV selectors optimize below the old h = 0.06 lower bound", {

  set.seed(123)
  samp <- rotasym::r_vMF(n = 500, mu = c(1, 0), kappa = 5000)
  h_rot <- bw_dir_rot(data = samp)
  expect_lt(h_rot, 0.06)
  expect_equal(bw_dir_lcv(data = samp)$h_opt, h_rot, tolerance = 0.5)
  expect_lt(bw_dir_lscv(data = samp)$h_opt, 0.06)

})

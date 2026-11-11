
## bw_dir_rot()

# The per-dimension closed forms of Proposition 2 in Garcia-Portugues (2013),
# as implemented before the single log-scale rewrite of bw_dir_rot().
# There is a typo in Garcia-Portugues (2013) in equation (6) for q = 2, where
# it displays an incorrect extra kappa:
# * WRONG:    kappa * (1 + 4 * kappa^2) * sinh(2 * kappa)
# * CORRECT:          (1 + 4 * kappa^2) * sinh(2 * kappa)
rot_naive <- function(data) {

  n <- nrow(data)
  q <- ncol(data) - 1
  kappa <- kappa_ml(data = data)
  if (q == 1) {

    num <- 4 * sqrt(pi) * besselI(nu = 0, x = kappa)^2
    den <- kappa * (2 * besselI(nu = 1, x = 2 * kappa) +
                      3 * kappa * besselI(nu = 2, x = 2 * kappa)) * n

  } else if (q == 2) {

    e <- exp(-2 * kappa)
    num <- 8 * (0.5 * (1 - e))^2
    den <- kappa * (-2 * kappa * (0.5 * (1 + e^2)) +
                      (1 + 4 * kappa^2) * (0.5 * (1 - e^2))) * n

  } else {

    num <- 4 * sqrt(pi) * besselI(nu = (q - 1) / 2, x = kappa,
                                  expon.scaled = TRUE)^2
    den <- kappa^((q + 1) / 2) *
      (2 * q * besselI(nu = (q + 1) / 2, x = 2 * kappa, expon.scaled = TRUE) +
         (2 + q) * kappa * besselI(nu = (q + 3) / 2, x = 2 * kappa,
                                   expon.scaled = TRUE)) * n

  }
  return((num / den)^(1 / (q + 4)))

}

test_that("bw_dir_rot reproduces the closed forms", {

  for (q in 1:4) {

    set.seed(200 + q)
    samp <- r_proj_norm(n = 80, mu = c(1, rep(0, q)), kappa = 4)
    expect_equal(bw_dir_rot(samp), rot_naive(samp), tolerance = 1e-8)

  }

})

test_that("bw_dir_rot is finite for large kappa and NA for degenerate data", {

  set.seed(5)
  for (true_kappa in c(1e5, 1e6)) {
    for (q in 1:3) {

      samp <- r_proj_norm(n = 50, mu = c(1, rep(0, q)), kappa = true_kappa)
      h <- bw_dir_rot(samp)
      expect_true(is.finite(h) && h > 0)

    }
  }

  # Identical observations have an infinite concentration
  ident <- matrix(c(0, 0, 1), nrow = 5, ncol = 3, byrow = TRUE)
  expect_warning(h <- bw_dir_rot(ident), "[Dd]egenerate")
  expect_true(is.na(h))

})

## bw_dir_lcv() and bw_dir_lscv()

# Von Mises--Fisher normalizing constant from base::besselI()
c_vmf_naive <- function(p, kappa) {

  return(kappa^(p / 2 - 1) /
           ((2 * pi)^(p / 2) * besselI(x = kappa, nu = p / 2 - 1)))

}

# Naive reference for the von Mises-Fisher LCV
lcv_naive <- function(data, h) {

  n <- nrow(data)
  p <- ncol(data)
  f <- sapply(h, function(hh) sapply(seq_len(n), function(i) {
    c_vmf_naive(p = p, kappa = 1 / hh^2) / (n - 1) *
      sum(exp(data[-i, ] %*% data[i, ] / hh^2))
  }))
  return(-colMeans(log(f)))

}

# Naive reference for the von Mises-Fisher LSCV
lscv_naive <- function(data, h) {

  n <- nrow(data)
  p <- ncol(data)
  return(sapply(h, function(hh) {

    ck <- c_vmf_naive(p = p, kappa = 1 / hh^2)
    ck2 <- c_vmf_naive(p = p, kappa = 2 / hh^2)
    s1 <- s2 <- 0
    for (i in 2:n) {
      for (j in 1:(i - 1)) {

        s1 <- s1 + exp(sum(data[i, ] * data[j, ]) / hh^2)
        s2 <- s2 + 1 / c_vmf_naive(p = p,
                                   kappa = norm2(data[i, ] + data[j, ]) /
                                     hh^2)

      }
    }
    -(4 * ck / (n * (n - 1)) * s1 - ck^2 / (n * ck2) - 2 * ck^2 / n^2 * s2)

  }))

}

test_that("The CV losses agree with their naive references", {

  h_grid <- exp(seq(log(0.4), log(1.2), l = 5))
  for (q in 1:3) {

    set.seed(7)
    samp <- r_proj_norm(n = 12, mu = c(1, rep(0, q)), kappa = 2)
    lcv <- suppressWarnings(bw_dir_lcv(data = samp, h_grid = h_grid,
                                       optim = FALSE))
    lscv <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                         optim = FALSE))
    expect_equal(lcv$CV_grid, lcv_naive(data = samp, h = h_grid),
                 tolerance = 1e-8)
    expect_equal(lscv$CV_grid, lscv_naive(data = samp, h = h_grid),
                 tolerance = 1e-7)

  }

})

test_that("The CV losses honor L and match vMF for L = exp(-z)", {

  # Likelihood CV
  set.seed(42)
  samp <- r_proj_norm(n = 30, mu = c(1, 0, 0), kappa = 2)
  h_grid <- exp(seq(log(0.1), log(1.2), l = 15))
  cv_f <- bw_dir_lcv(data = samp, h_grid = h_grid, optim = FALSE)
  cv_r <- bw_dir_lcv(data = samp, h_grid = h_grid, optim = FALSE,
                     L = function(z) exp(-z))
  expect_equal(cv_f$CV_grid, cv_r$CV_grid, tolerance = 1e-3)
  expect_equal(cv_f$h_opt, cv_r$h_opt)

  # Least squares CV
  set.seed(42)
  samp <- r_proj_norm(n = 20, mu = c(1, 0, 0), kappa = 2)
  h_grid <- c(0.4, 0.7, 1)
  vmf <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                      optim = FALSE))
  exp1 <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                       optim = FALSE, L = function(z) exp(-z)))
  exp2 <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                       optim = FALSE,
                                       L = function(z) exp(-2 * z)))
  expect_equal(exp1$CV_grid, vmf$CV_grid, tolerance = 1e-3)
  expect_false(isTRUE(all.equal(exp2$CV_grid, exp1$CV_grid)))

})

test_that("bw_dir_lscv with a custom kernel is finite for q = 1, 2, 3", {

  h_grid <- exp(seq(log(0.3), log(1.2), l = 6))
  for (q in 1:3) {

    set.seed(42)
    samp <- r_proj_norm(n = 20, mu = c(1, rep(0, q)), kappa = 2)
    lscv <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                         optim = FALSE,
                                         L = function(z) exp(-z)))
    expect_length(lscv$CV_grid, length(h_grid))
    expect_true(all(is.finite(lscv$CV_grid)))

  }

})

test_that("The CV losses accumulate the pairwise sums over several blocks", {

  # n = 1001 splits bw_dir_lscv() into blocks, and n = 5001 > 5000 also stops
  # bw_dir_lcv() from caching the whole matrix of inner products.
  set.seed(7)
  samp <- r_proj_norm(n = 1001, mu = c(1, 0, 0), kappa = 4)
  h_grid <- c(0.1, 0.3, 0.8)
  lscv <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                       optim = FALSE))
  lscv_r <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                         optim = FALSE,
                                         L = function(z) exp(-z)))
  expect_equal(lscv$CV_grid, lscv_r$CV_grid, tolerance = 1e-4)
  set.seed(7)
  samp <- r_proj_norm(n = 5001, mu = c(1, 0, 0), kappa = 4)
  h_grid <- c(0.2, 0.5, 1)
  lcv <- suppressWarnings(bw_dir_lcv(data = samp, h_grid = h_grid,
                                     optim = FALSE))
  lcv_r <- suppressWarnings(bw_dir_lcv(data = samp, h_grid = h_grid,
                                       optim = FALSE, L = function(z) exp(-z)))
  expect_equal(lcv$CV_grid, lcv_r$CV_grid, tolerance = 1e-4)

})

test_that("The CV selectors optimize consistently, also below h = 0.06", {

  # optim() and the grid search agree
  set.seed(42)
  samp <- r_proj_norm(n = 40, mu = c(1, 0, 0), kappa = 2)
  expect_equal(bw_dir_lcv(data = samp, optim = TRUE)$h_opt,
               bw_dir_lcv(data = samp, optim = FALSE)$h_opt, tolerance = 0.1)
  expect_equal(bw_dir_lscv(data = samp, optim = TRUE)$h_opt,
               bw_dir_lscv(data = samp, optim = FALSE)$h_opt, tolerance = 0.1)

  # Concentrated data, with bandwidths below the old h = 0.06 lower bound
  set.seed(123)
  samp <- r_proj_norm(n = 500, mu = c(1, 0), kappa = 5000)
  h_rot <- bw_dir_rot(data = samp)
  expect_lt(h_rot, 0.06)
  expect_equal(bw_dir_lcv(data = samp)$h_opt, h_rot, tolerance = 0.5)
  expect_lt(bw_dir_lscv(data = samp)$h_opt, 0.06)

})

test_that("The CV losses are finite for n = 2, antipodal pairs, and tiny h", {

  # Two observations leave a single leave-one-out term; an antipodal pair gives
  # a dangerous ||X_i + X_j|| = 0, which old code divided by.
  samp2 <- rbind(c(1, 0, 0), c(0, 1, 0))
  anti <- rbind(c(1, 0, 0), c(-1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  h_grid <- c(0.5, 1)
  expect_true(all(is.finite(suppressWarnings(
    bw_dir_lcv(data = samp2, h_grid = h_grid, optim = FALSE)$CV_grid))))
  expect_true(all(is.finite(suppressWarnings(
    bw_dir_lscv(data = samp2, h_grid = h_grid, optim = FALSE)$CV_grid))))
  expect_true(all(is.finite(suppressWarnings(
    bw_dir_lscv(data = anti, h_grid = h_grid, optim = FALSE)$CV_grid))))

  # Bandwidths down to sqrt(.Machine$double.eps)
  h_grid <- c(sqrt(.Machine$double.eps), 1e-6, 1e-3, 0.01, 0.05, 0.5)
  for (q in 1:3) {

    set.seed(7)
    samp <- r_proj_norm(n = 20, mu = c(1, rep(0, q)), kappa = 3)
    lcv <- suppressWarnings(bw_dir_lcv(data = samp, h_grid = h_grid,
                                       optim = FALSE)$CV_grid)
    lscv <- suppressWarnings(bw_dir_lscv(data = samp, h_grid = h_grid,
                                         optim = FALSE)$CV_grid)
    expect_true(all(is.finite(lcv)))
    expect_true(all(is.finite(lscv)))

  }

})

## bic_vmf_mix()

test_that("bic_vmf_mix returns the documented structure", {

  set.seed(42)
  samp <- r_proj_norm(n = 60, mu = c(1, 0, 0), kappa = 3)
  out <- capture.output(fit <- bic_vmf_mix(data = samp, verbose = TRUE,
                                           plot_it = TRUE))
  expect_true(any(grepl("Fitting mixtures", out)))
  expect_named(fit, c("best_fit", "fits", "BICs"))
  expect_named(fit$best_fit, c("mu_hat", "kappa_hat", "p_hat"))
  expect_equal(sum(fit$best_fit$p_hat), 1, tolerance = 1e-6)

})

test_that("bic_vmf_mix accepts crit as a function", {

  set.seed(42)
  samp <- r_proj_norm(n = 60, mu = c(1, 0, 0), kappa = 3)
  set.seed(1)
  fit_name <- bic_vmf_mix(data = samp, crit = "AIC", iterative = FALSE)
  set.seed(1)
  fit_fun <- bic_vmf_mix(data = samp, crit = AIC, iterative = FALSE)
  expect_equal(fit_name$BICs, fit_fun$BICs)
  set.seed(1)
  AICc <- function(object) {

    ll <- logLik(object)
    k <- attr(ll, "df")
    n <- attr(ll, "nobs")
    return(AIC(object) + 2 * k * (k + 1) / (n - k - 1))

  }
  fit_aicc <- bic_vmf_mix(data = samp, crit = AICc, iterative = FALSE)
  k <- sapply(fit_name$fits,
              function(x) ifelse(is.na(x[1]), NA, attr(logLik(x), "df")))
  expect_equal(fit_aicc$BICs,
               fit_name$BICs + 2 * k * (k + 1) / (nrow(samp) - k - 1),
               ignore_attr = TRUE)

})

test_that("bic_vmf_mix runs the iterative extension of the search", {

  set.seed(42)
  q <- 2
  samp <- rbind(
    r_proj_norm(n = 30, mu = c(rep(0, q), 1), kappa = 8),
    r_proj_norm(n = 30, mu = c(rep(0, q), -1), kappa = 8),
    r_proj_norm(n = 30, mu = c(1, rep(0, q)), kappa = 8))
  fit <- bic_vmf_mix(data = samp, M_bound = 2, M_neig = 3, iterative = TRUE)
  expect_named(fit, c("best_fit", "fits", "BICs"))
  expect_true(length(fit$BICs) > 2)
  expect_gte(length(fit$best_fit$kappa_hat), 1)

})

## Resistance to degenerate mixtures

test_that("kappa_max screens collapsed components only", {

  # A component collapsing onto a few of the thirty observations reaches a
  # concentration above 1e4 and gets selected.
  set.seed(57)
  samp <- r_proj_norm(n = 30, mu = c(1, 0), kappa = 3)
  set.seed(1)
  fit <- bic_vmf_mix(data = samp)
  expect_true(all(fit$best_fit$kappa_hat <= 1e4))
  # bw_dir_ami inherits the component names of movMF's alpha, hence ignore_attr
  expect_equal(bw_dir_ami(data = samp, fit_mix = fit), bw_dir_rot(data = samp),
               tolerance = 1e-4, ignore_attr = TRUE)

  # kappa_max = Inf disables the screening, recovering the degenerate mixture
  set.seed(1)
  fit_inf <- bic_vmf_mix(data = samp, kappa_max = Inf)
  expect_true(any(fit_inf$best_fit$kappa_hat > 1e4))

  # A genuinely concentrated sample is kept: its single-component fit has kappa
  # about 5000, well below kappa_max.
  set.seed(123)
  samp <- r_proj_norm(n = 500, mu = c(1, 0), kappa = 5000)
  set.seed(1)
  fit <- suppressWarnings(bic_vmf_mix(data = samp))
  expect_true(all(is.finite(fit$best_fit$kappa_hat)))
  expect_gt(max(fit$best_fit$kappa_hat), 1e3)

  # So is a multimodal one, whose pooled resultant length is nearly zero
  set.seed(1)
  samp <- rbind(r_proj_norm(n = 50, mu = c(1, 0, 0), kappa = 20),
                r_proj_norm(n = 50, mu = c(-1, 0, 0), kappa = 20))
  set.seed(1)
  fit <- bic_vmf_mix(data = samp)
  expect_length(fit$best_fit$p_hat, 2)
  expect_true(all(fit$best_fit$kappa_hat > 5))

})

test_that("bic_vmf_mix falls back to a single vMF and discards worse refits", {

  # All the movMF() fits failing leaves no criterion to minimize. The fallback
  # returns mu_hat as a one-row matrix, as the plug-in selectors expect, also
  # when plotting the criterion.
  set.seed(1)
  samp <- r_proj_norm(n = 20, mu = c(1, 0, 0), kappa = 2)
  fail <- function(...) stop("no fit")
  expect_warning(fit <- with_mocked_bindings(code = bic_vmf_mix(data = samp,
                                                                plot_it = TRUE),
                                             movMF = fail, .package = "movMF"),
                 "single vMF")
  expect_equal(fit$best_fit$p_hat, 1)
  expect_equal(dim(fit$best_fit$mu_hat), c(1L, 3L))
  expect_equal(fit$best_fit$kappa_hat, kappa_ml(data = samp))
  expect_true(all(is.na(fit$BICs)))

  # A single vMF makes bw_dir_ami() coincide with the rule of thumb
  expect_equal(bw_dir_ami(data = samp, fit_mix = fit),
               bw_dir_rot(data = samp))

  # The last movMF() call refines the selected mixture. Forcing it to return a
  # single component makes it worse than the selected two, so it is discarded.
  set.seed(1)
  samp <- rbind(r_proj_norm(n = 15, mu = c(1, 0, 0), kappa = 10),
                r_proj_norm(n = 15, mu = c(-1, 0, 0), kappa = 10))
  M_bound <- ceiling(log(nrow(samp)))
  real_movMF <- movMF::movMF
  calls <- 0
  worse <- function(x, k, ...) {

    calls <<- calls + 1
    return(real_movMF(x = x, k = ifelse(calls > M_bound, 1, k), ...))

  }
  set.seed(1)
  fit <- with_mocked_bindings(code = bic_vmf_mix(data = samp,
                                                 iterative = FALSE),
                              movMF = worse, .package = "movMF")
  expect_length(fit$best_fit$p_hat, 2)

})

test_that("The plug-in selectors survive extreme mixture concentrations", {

  # The curvature is finite and increasing in the concentration up to 1e8
  for (q in 1:3) {

    mu <- rbind(c(rep(0, q), 1))
    R <- sapply(c(1e4, 1e5, 1e6, 1e8), function(k)
      R_Psi_mixvmf(q = q, mu = mu, kappa = k, p = 1))
    expect_true(all(is.finite(R)) && all(R > 0))
    expect_true(all(diff(R) > 0))

  }

  # A mixture with a component of concentration 1.2e5
  mu <- rbind(c(1, 0, 0), c(0, 1, 0))
  kappa <- c(4, 1.2e5)
  p <- c(0.9, 0.1)
  R <- R_Psi_mixvmf(q = 2, mu = mu, kappa = kappa, p = p)
  expect_true(is.finite(R) && R > 0)
  set.seed(42)
  samp <- r_proj_norm(n = 20, mu = c(1, 0, 0), kappa = 4)
  fit_mix <- list(best_fit = list(mu_hat = mu, kappa_hat = kappa, p_hat = p))
  h_ami <- bw_dir_ami(data = samp, fit_mix = fit_mix)
  expect_true(is.finite(h_ami) && h_ami > 0)
  set.seed(42)
  mise <- suppressWarnings(bw_dir_emi(data = samp, fit_mix = fit_mix,
                                      optim = FALSE, plot_it = FALSE,
                                      h_grid = c(0.2, 0.5, 1))$MISE_grid)
  expect_true(all(is.finite(mise)))

  # A collapsed component of concentration 1e18, capped in the proposal
  set.seed(987202226)
  q <- 6
  data <- normalize(matrix(rnorm(n = 300 * (q + 1)), nrow = 300, ncol = q + 1))
  mu <- normalize(matrix(rnorm(n = 3 * (q + 1)), nrow = 3, ncol = q + 1))
  fit_mix <- list(best_fit = list(mu_hat = mu, kappa_hat = c(50, 500, 1e18),
                                  p_hat = c(0.4, 0.4, 0.2)))
  expect_warning(emi <- bw_dir_emi(data = data, fit_mix = fit_mix,
                                   optim = TRUE, plot_it = FALSE),
                 "collapsed")
  expect_true(is.finite(emi$h_opt) && emi$h_opt > 0)

})

## R_Psi_mixvmf()

test_that("R_Psi_mixvmf numeric route matches the exact one", {

  q <- c(1, 2, 3, 5)
  kappa <- c(3, 4, 3, 500)
  tol <- c(1e-2, 1e-2, 5e-2, 5e-2)
  for (i in seq_along(q)) {

    mu <- c(1, rep(0, q[i]))
    set.seed(987202226)
    R_exact <- R_Psi_mixvmf(q = q[i], mu = rbind(mu), kappa = kappa[i], p = 1)
    R_numer <- R_Psi_mixvmf(q = q[i], mu = rbind(mu, mu),
                            kappa = rep(kappa[i], 2), p = c(0.5, 0.5))
    expect_equal(R_numer, R_exact, tolerance = tol[i])

  }

})

test_that("R_Psi_mixvmf checks argument dimensions", {

  expect_error(R_Psi_mixvmf(q = 2, mu = rbind(c(0, 0, 1)), kappa = c(1, 2),
                            p = c(0.5, 0.5)), "dimension")

})

## bw_dir_ami()

test_that("bw_dir_ami returns a positive bandwidth on S^1, S^2, and S^3", {

  for (q in 1:3) {

    set.seed(42)
    samp <- r_proj_norm(n = 40, mu = c(1, rep(0, q)), kappa = 3)
    ami <- bw_dir_ami(data = samp)
    expect_true(is.finite(ami) && ami > 0)

  }

})

## bw_dir_emi()

test_that("bw_dir_emi returns a positive bandwidth on S^1, S^2, and S^3", {

  for (q in 1:3) {

    set.seed(42)
    samp <- r_proj_norm(n = 40, mu = c(1, rep(0, q)), kappa = 3)
    emi <- bw_dir_emi(data = samp, optim = TRUE)
    expect_true(is.finite(emi$h_opt) && emi$h_opt > 0)

  }

})

test_that("bw_dir_emi MISE is finite on a grid extending below h = 0.05", {

  h_grid <- exp(seq(log(0.03), log(1.5), l = 15))
  for (q in 1:3) {

    set.seed(7)
    samp <- r_proj_norm(n = 30, mu = c(1, rep(0, q)), kappa = 3)
    fit_mix <- list(best_fit = list(mu_hat = rbind(mu_ml(samp)),
                                    kappa_hat = 3, p_hat = 1))
    set.seed(7)
    emi <- suppressWarnings(bw_dir_emi(samp, fit_mix = fit_mix, optim = FALSE,
                                       plot_it = TRUE, h_grid = h_grid))
    expect_length(emi$MISE_grid, length(h_grid))
    expect_true(all(is.finite(emi$MISE_grid)) && all(emi$MISE_grid > 0))
    expect_true(is.finite(emi$h_opt))

  }

})

# A mixture concentrated enough on S^5 for uniform draws to miss its mass
emi_case <- function() {

  set.seed(987202226)
  q <- 5
  mu <- normalize(matrix(rnorm(n = 2 * (q + 1)), nrow = 2, ncol = q + 1))
  kappa <- c(300, 800)
  data <- rbind(r_proj_norm(n = 50, mu = mu[1, ], kappa = kappa[1]),
                r_proj_norm(n = 50, mu = mu[2, ], kappa = kappa[2]))
  return(list(data = data, fit_mix = list(best_fit =
    list(mu_hat = mu, kappa_hat = kappa, p_hat = c(0.5, 0.5)))))

}

test_that("bw_dir_emi has an interior optimum if concentrated", {

  cs <- emi_case()
  h_grid <- exp(seq(log(0.01), log(1), l = 12))
  emi <- suppressWarnings(bw_dir_emi(data = cs$data, fit_mix = cs$fit_mix,
                                     optim = FALSE, plot_it = FALSE,
                                     h_grid = h_grid))
  expect_true(all(is.finite(emi$MISE_grid)))
  ind_min <- which.min(emi$MISE_grid)
  expect_gt(ind_min, 1)
  expect_lt(ind_min, length(h_grid))
  expect_lt(emi$MISE_grid[ind_min], 0.9 * emi$MISE_grid[length(h_grid)])

})

test_that("bw_dir_emi draws its sample once per call, so optim() is stable", {

  cs <- emi_case()
  h_opt <- sapply(1:3, function(s) {
    set.seed(s)
    suppressWarnings(bw_dir_emi(data = cs$data, fit_mix = cs$fit_mix,
                                optim = TRUE)$h_opt)
  })
  expect_true(all(is.finite(h_opt) & h_opt > 0))
  expect_lt(diff(range(h_opt)) / mean(h_opt), 0.1)

})

## Grid boundary warnings

test_that("The selectors warn when the optimum is at the extreme of h_grid", {

  set.seed(42)
  samp <- r_proj_norm(n = 20, mu = c(1, 0, 0), kappa = 2)
  fit_mix <- list(best_fit = list(mu_hat = rbind(c(1, 0, 0)), kappa_hat = 2,
                                  p_hat = 1))
  expect_warning(bw_dir_lcv(data = samp, h_grid = c(1.4, 1.5), optim = FALSE),
                 "extreme")
  expect_warning(bw_dir_lscv(data = samp, h_grid = c(1.4, 1.5), optim = FALSE),
                 "extreme")
  expect_warning(bw_dir_emi(data = samp, fit_mix = fit_mix, optim = FALSE,
                            h_grid = c(1.4, 1.5)), "extreme")

})

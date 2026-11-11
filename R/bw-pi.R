

#' @title Fitting mixtures of von Mises--Fisher distributions
#'
#' @description Fitting mixtures of von Mises--Fisher distributions by the
#' Expectation-maximization algorithm, with determination of the optimal number
#' of mixture components.
#'
#' @inheritParams kde_dir
#' @param M_bound bound for the number of components in the mixtures. If it is
#' not enough, the search for the mixture with minimum \code{crit} will continue
#' from \code{M_bound + 1} if \code{iterative = TRUE}. Defaults to
#' \code{ceiling(log(nrow(data)))}.
#' @param M_neig number of neighbors explored around the optimal number of
#' mixture components. Defaults to \code{3}.
#' @param crit information criterion to be minimized, given as the name of a
#' function or as a function itself. It must take a fitted mixture and return a
#' single number. Defaults to \code{"BIC"}; see the examples for a custom
#' criterion.
#' @param iterative keep exploring higher number of components if the optimum is
#' attained at \code{M_bound}? Defaults to \code{TRUE}.
#' @inheritParams bw_dir_cv
#' @param verbose display fitting progress? Defaults to \code{FALSE}.
#' @param kappa_max maximum concentration allowed in a fitted mixture. Mixtures
#' with a component above it are discarded. Set to \code{Inf} to disable the
#' screening. Defaults to \code{1e4}.
#' @return A list with entries:
#' \itemize{
#'   \item \code{best_fit}: a list with estimated mixture parameters
#'   \code{mu_hat}, \code{kappa_hat}, and \code{p_hat} of the best-fitting
#'   mixture according to \code{crit}.
#'   \item \code{fits}: a list with the fitted mixtures.
#'   \item \code{BICs}: a vector with the BICs (or other information criterion)
#'   of the fitted mixtures.
#' }
#' @details
#' See Algorithm 3 in García-Portugués (2013). The Expectation-Maximization
#' fit is performed with \code{\link[movMF]{movMF}}.
#'
#' The \code{kappa_max} screening allows discarding fits with spurious
#' components featuring high concentration.
#' @references
#' García-Portugués, E. (2013). Exact risk improvement of bandwidth selectors
#' for kernel density estimation with directional data. \emph{Electronic
#' Journal of Statistics}, 7:1655--1685.
#' \doi{10.1214/13-EJS821}
#'
#' Hornik, K. and Grün, B. (2014). movMF: An R Package for Fitting Mixtures of
#' von Mises--Fisher Distributions. \emph{Journal of Statistical Software},
#' 58(10):1--31. \doi{10.18637/jss.v058.i10}
#' @seealso \code{\link{bw_dir_pi}} for the plug-in bandwidth selectors that
#' build on this mixture fit, and \code{\link{vmf}} for the single-component
#' estimators used as a fallback.
#' @examples
#' # Custom criterion: the AIC with the finite-sample correction
#' # 2 * k * (k + 1) / (n - k - 1), for k parameters and n observations.
#' AICc <- function(object) {
#'   ll <- logLik(object)
#'   k <- attr(ll, "df")
#'   n <- attr(ll, "nobs")
#'   AIC(object) + ifelse(n - k - 1 > 0, 2 * k * (k + 1) / (n - k - 1), Inf)
#' }
#'
#' # Sample from a mixture of three von Mises--Fisher distributions
#' q <- 2
#' n <- 300
#' if (requireNamespace("rotasym", quietly = TRUE)) {
#'
#'   set.seed(42)
#'   samp <- rbind(rotasym::r_vMF(n = n / 3, mu = c(rep(0, q), 1), kappa = 5),
#'                 rotasym::r_vMF(n = n / 3, mu = c(rep(0, q), -1), kappa = 5),
#'                 rotasym::r_vMF(n = n / 3, mu = c(1, rep(0, q)), kappa = 5))
#'
#'   # Mixture fit
#'   print(bic_vmf_mix(data = samp, plot_it = TRUE, verbose = TRUE))
#'   print(bic_vmf_mix(data = samp, crit = AICc, plot_it = TRUE)$BICs)
#'
#' }
#' @export
bic_vmf_mix <- function(data, M_bound = ceiling(log(nrow(data))), M_neig = 3,
                        crit = "BIC", iterative = TRUE, plot_it = FALSE,
                        verbose = FALSE, kappa_max = 1e4) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Create list for mixtures
  fit_mixs <- vector("list", length = M_bound)

  # Create vector for BICs
  BICs <- rep(NA, M_bound)

  # For computing the information criterion, given as a function name or as a
  # function itself. The name is only kept for labeling the plot.
  crit_name <- if (is.character(crit)) crit else deparse(substitute(crit))[1]
  crit <- match.fun(crit)

  # Mixture fit with k components, NA if movMF() fails
  fit_k <- function(k, maxiter = 100) {

    tryCatch(movMF::movMF(x = data, k = k, maxiter = maxiter, start = "S"),
             error = function(e) NA)

  }

  # Information criterion of a fit, NA if the fit failed or if it is degenerate
  # (a component above kappa_max).
  get_BIC <- function(x) {

    if (is.na(x[1]) || any(norm2(x$theta) > kappa_max)) NA else crit(x)

  }

  # Mixture fitting for 1:M_bound
  if (verbose) {

    cat("Fitting mixtures...\n")
    pb <- txtProgressBar(style = 3)

  }
  for (k in 1:M_bound) {

    # Fit mixture
    fit_mixs[[k]] <- fit_k(k = k)
    BICs[k] <- get_BIC(fit_mixs[[k]])
    if (verbose) setTxtProgressBar(pb, k / M_bound)

  }

  # Seek for the minimum. Note that NA and NaN values are discarded
  M_opt <- which.min(BICs)

  # What if M_opt is near the maximum of the grid? Run the iterative procedure.
  # M_opt is empty if no mixture could be fitted, in which case the single von
  # Mises--Fisher fallback below takes over.
  if (iterative && length(M_opt) == 1) {

    if (verbose) cat("\nComputing iterative procedure...")

    # Check if M_opt is in (M_bound - M_neig + 1):M_bound
    new_M_bound <- M_bound
    if ((M_bound - M_neig + 1) <= M_opt) {

      # Increase new_M_bound before loop
      new_M_bound <- M_bound + 1

      # Loop until ensuring the minimum is within the new extended range
      repeat {

        # New estimation, added to the previous ones
        fit_mixs_new <- fit_k(k = new_M_bound)
        fit_mixs <- c(fit_mixs, list(fit_mixs_new))
        BICs <- c(BICs, get_BIC(fit_mixs_new))

        # Seek for the minimum. Note that NA and NaN values are discarded
        M_opt <- which.min(BICs)

        # Stop iteration if M_opt < new_M_bound - M_neig + 1; otherwise
        # increase the bound.
        if (M_opt < (new_M_bound - M_neig + 1)) {

          break

        } else {

          new_M_bound <- new_M_bound + 1

        }

      }

      # Set the old bound to the obtained new_M_bound
      M_bound <- new_M_bound

    }

    if (verbose) cat("Done.\n")

  }

  # Extract parameters
  if (length(M_opt) == 0) {

    # Maximum likelihood fit of a vMF if everything fails
    warning("No mixture could be fitted; returning a single vMF fit.")
    mu_hat <- rbind(mu_ml(data = data))
    kappa_hat <- kappa_ml(data = data)
    p_hat <- 1

  } else {

    # Compute a more precise solution for the optimum, taken only if there are
    # no convergence issues, it does not worsen crit, and it is not degenerate,
    # as a longer run may land on a spiked mixture.
    fit_mix <- fit_mixs[[M_opt]]
    fit_mixs_better <- fit_k(k = M_opt, maxiter = 300)
    BIC_better <- get_BIC(fit_mixs_better)
    if (!is.na(BIC_better) && BIC_better <= BICs[M_opt]) {

      fit_mix <- fit_mixs_better

    }

    # Extract parameters
    kappa_hat <- norm2(fit_mix$theta)
    mu_hat <- fit_mix$theta / kappa_hat
    p_hat <- fit_mix$alpha

  }

  # Display a graph? Not if no mixture could be fitted, as there is nothing to
  # plot.
  if (plot_it && any(is.finite(BICs))) {

    plot(seq_along(BICs), BICs, xlab = "M", ylab = crit_name,
         main = paste("Min", crit_name), pch = 19, cex = 1.5)
    points(M_opt, BICs[M_opt], pch = 19, cex = 1.5, col = 2)

  }

  # Final object
  return(list(best_fit =
                list(mu_hat = mu_hat, kappa_hat = kappa_hat, p_hat = p_hat),
              fits = fit_mixs, BICs = BICs))

}


#' @title Plug-in bandwidth selectors for directional data
#'
#' @description Plug-in bandwidth selectors for kernel density estimation with
#' directional data, including Rule-Of-Thumb (ROT), Asymptotic MIxtures (AMI),
#' and Exact MIxtures (EMI).
#'
#' @inheritParams kde_dir
#' @inheritParams bic_vmf_mix
#' @inheritParams bw_dir_cv
#' @inheritParams vmf
#' @param fit_mix output from \code{\link{bic_vmf_mix}}. Computed internally if
#' \code{NULL} (default).
#' @return
#' \itemize{
#'   \item \code{bw_dir_rot}, \code{bw_dir_ami}: the selected bandwidth, a
#'   scalar.
#'   \item \code{bw_dir_emi}: a list with entries \code{h_opt} (selected
#'   bandwidth), \code{h_grid} (\code{h_grid}, if used, otherwise
#'   \code{NULL}), \code{MISE_opt} (minimum of the MISE loss), and
#'   \code{MISE_grid} (value of the MISE function at \code{h_grid}, if used,
#'   otherwise \code{NULL}).
#'   \item \code{R_Psi_mixvmf}: the curvature term of the mixture, a scalar.
#' }
#' @details
#' See Algorithms 1 (AMI) and 2 (EMI) in García-Portugués (2013). The ROT
#' selector is implemented according to Proposition 2, \bold{but} without
#' the paper's typo in equation (6), case \eqn{q = 2}, where an incorrect
#' extra \eqn{\hat\kappa} appears premultiplying
#' \eqn{(1 + 4 \hat\kappa^2) \sinh(2 \hat\kappa)} in the denominator.
#'
#' \code{bw_dir_ami} uses \code{R_Psi_mixvmf} for computing the curvature
#' term of a mixture of von Mises--Fisher densities.
#'
#' \code{bw_dir_emi} employs Monte Carlo integration for \eqn{q > 2}, which
#' results in a random output. Use \code{set.seed} before to avoid it.
#'
#' The selectors warn and return \code{NA} (\code{h_opt = NA} for
#' \code{bw_dir_emi}) for degenerate concentrations.
#' @section Warning:
#' \code{bw_dir_emi} overflows for large dimensions \eqn{q \approx 250}, and
#' earlier for concentrated data.
#' @references
#' García-Portugués, E. (2013). Exact risk improvement of bandwidth selectors
#' for kernel density estimation with directional data. \emph{Electronic
#' Journal of Statistics}, 7:1655--1685.
#' \doi{10.1214/13-EJS821}
#' @seealso \code{\link{kde_dir}} for the estimator, \code{\link{bw_dir_cv}} for
#' cross-validation alternatives, and \code{\link{bic_vmf_mix}} for the mixture
#' fit these selectors rely on.
#' @examples
#' # Sample
#' n <- 25
#' q <- 2
#' if (requireNamespace("rotasym", quietly = TRUE)) {
#'
#'   set.seed(42)
#'   samp <- rotasym::r_vMF(n = n, mu = c(1, rep(0, q)), kappa = 2)
#'
#'   # Mixture fit
#'   fit_mix <- bic_vmf_mix(data = samp, plot_it = TRUE)
#'
#'   # ROT
#'   print(bw_dir_rot(samp))
#'
#'   # AMI
#'   print(bw_dir_ami(samp))
#'   print(bw_dir_ami(samp, fit_mix = fit_mix))
#'   print(bw_dir_ami(samp, fit_mix = fit_mix, L = function(x) exp(-x)))
#'
#'   # EMI
#'   print(bw_dir_emi(samp))
#'   print(bw_dir_emi(samp, fit_mix = fit_mix, optim = FALSE, plot_it = TRUE))
#'
#' }
#' @name bw_dir_pi


#' @rdname bw_dir_pi
#' @export
bw_dir_rot <- function(data) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Rule of thumb of Proposition 2 in Garcia-Portugues (2013): the AMI selector
  # with a single von Mises--Fisher as reference, whose curvature
  # R_Psi_mixvmf() evaluates in log scale for arbitrarily large concentrations.
  # There is a typo in Garcia-Portugues (2013) in equation (6) for q = 2, where
  # it displays an incorrect extra kappa:
  # * WRONG:    kappa * (1 + 4 * kappa^2) * sinh(2 * kappa)
  # * CORRECT:          (1 + 4 * kappa^2) * sinh(2 * kappa)
  # This is not relevant for the general case q >= 1 computed here.
  fit_mix <- list(best_fit = list(mu_hat = rbind(mu_ml(data = data)),
                                  kappa_hat = kappa_ml(data = data),
                                  p_hat = 1))

  # Result
  return(bw_dir_ami(data = data, fit_mix = fit_mix))

}


#' @rdname bw_dir_pi
#' @export
bw_dir_ami <- function(data, fit_mix = NULL, L = NULL) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Constants needed for h_AMISE
  n <- nrow(data)
  q <- ncol(data) - 1
  lambda <- lambda_L(L = L, q = q)
  b <- b_L(L = L, q = q)
  d <- d_L(L = L, q = q)

  # Fit mixtures
  if (is.null(fit_mix)) {

    fit_mix <- bic_vmf_mix(data = data)

  }
  mu <- fit_mix$best_fit$mu_hat
  kappa <- fit_mix$best_fit$kappa_hat
  p <- fit_mix$best_fit$p_hat

  # Degenerate concentrations (infinite or all zero)
  if (any(!is.finite(kappa)) || all(kappa <= 0)) {

    warning(paste0("Degenerate concentration estimates (kappa = ",
                   paste(kappa, collapse = ", "), "); returning NA."))
    return(NA)

  }

  # Computation of curvature
  R <- R_Psi_mixvmf(q = q, mu = mu, kappa = kappa, p = p)

  # h_AMISE
  h <- ((q * d) / (4 * b^2 * lambda * R * n))^(1 / (q + 4))
  return(h)

}


#' @rdname bw_dir_pi
#' @export
R_Psi_mixvmf <- function(q, mu, kappa, p) {

  # Dimension checks
  if (length(p) != nrow(mu) || length(p) != length(kappa)) {

    stop("Check the dimension of the arguments")

  }

  # Psi operator for a single vMF(mu, kappa)
  Psi_vmf <- function(x, mu, kappa) {

    x_mu <- drop(x %*% mu)
    log_const <- log_c_vmf(q = q, kappa = kappa)
    return(kappa * exp(log_const + kappa * x_mu) *
             (-x_mu + kappa / q * (1 - x_mu^2)))

  }

  # Exact expression
  if (length(p) == 1) {

    # The uniform density has zero curvature
    if (kappa == 0) {

      R <- 0

    } else {

      # Evaluated in log scale (with a log-sum-exp for the two-term sum), so
      # that the curvature stays finite for large concentrations. The
      # exp(2 * kappa) factors of const and R cancel.
      log_const <- ((q + 1) / 2) * log(kappa) -
        ((q + 2) * log(2) + ((q + 1) / 2) * log(pi) + log(q) +
           2 * log_besselI_scaled(nu = (q - 1) / 2, x = kappa))
      log_1 <- log(2 * q) +
        log_besselI_scaled(nu = (q + 1) / 2, x = 2 * kappa)
      log_2 <- log((2 + q) * kappa) +
        log_besselI_scaled(nu = (q + 3) / 2, x = 2 * kappa)
      R <- exp(log_const + log_sum_exp(c(log_1, log_2)))

    }

  # Numeric approximation
  } else {

    # Curvature integrand
    integrand <- function(x) {

      rowSums(sapply(seq_along(p), function(i)
        Psi_vmf(x = x, mu = mu[i, ], kappa = kappa[i])
      ) %*% p)^2

    }

    # Integration, by importance sampling for q > 2. The proposal is a
    # defensive mixture of the mixture itself with a uniform, which keeps the
    # importance weights bounded where the mixture is close to zero.
    prop <- switch((q > 2) + 1, NULL,
                   safe_mixvmf(q = q, mu = mu, kappa = kappa, p = p, M = 1e4))
    R <- int_dir(f = integrand, q = q, sim = prop$sim, f_sim = prop$f_sim)

  }

  # Curvature term
  return(R)

}


#' @rdname bw_dir_pi
#' @export
bw_dir_emi <- function(data, fit_mix = NULL, optim = TRUE,
                       h_grid = exp(seq(log(0.05), log(1.5), l = 100)),
                       plot_it = FALSE, optim_par = 0.25,
                       optim_lower = max(1e-4, exp(-700 / (ncol(data) - 1)) /
                                           sqrt(2 * pi)),
                       optim_upper = 10) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Dimensions
  n <- nrow(data)
  q <- ncol(data) - 1

  # Fit mixtures
  if (is.null(fit_mix)) {

    fit_mix <- bic_vmf_mix(data = data)

  }
  mu <- fit_mix$best_fit$mu_hat
  kappa <- fit_mix$best_fit$kappa_hat
  p <- fit_mix$best_fit$p_hat

  # Infinite concentrations make the MISE undefined
  if (any(!is.finite(kappa))) {

    warning(paste0("Degenerate concentration estimates (kappa = ",
                   paste(kappa, collapse = ", "), "); returning NA."))
    return(list(h_opt = NA, h_grid = NULL, MISE_opt = NA, MISE_grid = NULL))

  }

  # Integrand for MISE function
  integrand <- function(x, h) {

    # Expectation of the estimator, the mixture over j of
    # c_q(kappa_j) * c_q(a) / c_q(r_j), with a = 1 / h^2 and
    # r_j = ||a * x + kappa_j * mu_j||. In terms of the exponentially scaled
    # constants \tilde c_q(kappa) = exp(kappa) * c_q(kappa), the remaining
    # exponent r_j - a - kappa_j is computed as
    #   r_j - a - kappa_j = (r_j^2 - (a + kappa_j)^2) / (r_j + a + kappa_j)
    #     = -2 * a * kappa_j * (1 - x'mu_j) / (r_j + a + kappa_j)
    a <- 1 / h^2
    log_c_a <- log_c_vmf(q = q, kappa = a, scaled = TRUE)
    log_Ef <- matrix(sapply(seq_along(p), function(j) {

      t_j <- drop(x %*% mu[j, ])
      r_j <- norm2(sweep(a * x, 2, kappa[j] * mu[j, ], "+"))
      log(p[j]) + log_c_vmf(q = q, kappa = kappa[j], scaled = TRUE) + log_c_a -
        log_c_vmf(q = q, kappa = r_j, scaled = TRUE) -
        2 * a * kappa[j] * (1 - t_j) / (r_j + a + kappa[j])

    }), nrow = nrow(x))
    Ef <- rowSums(exp(log_Ef))

    # Mixture of von Mises--Fisher
    f_mix <- d_mixvmf(x = x, mu = mu, kappa = kappa, p = p)

    # Squared bias
    Bias2 <- (Ef - f_mix)^2

    # Second term of variance, the first term can be exactly integrated
    Var1 <- -Ef^2 / n

    # Integrand without the first variance term
    res <- Bias2 + Var1
    return(res)

  }

  # Drawing the sample once makes mise_star() deterministic in h. For q > 2, the
  # proposal is a defensive mixture of the fitted mixture with a uniform, as
  # in R_Psi_mixvmf().
  prop <- switch((q > 2) + 1, NULL,
                 safe_mixvmf(q = q, mu = mu, kappa = kappa, p = p, M = 1e4))

  # Definition of mise_star: the integral plus the exactly integrated first
  # variance term, c_h(h)^2 / (n * c_h(h / sqrt(2))) as in bw_dir_lscv().
  mise_star <- function(h) {

    int <- int_dir(f = integrand, q = q, sim = prop$sim, f_sim = prop$f_sim,
                   h = h)
    return(int + exp(2 * c_h(h = h, q = q, log = TRUE) -
                       c_h(h = h / sqrt(2), q = q, log = TRUE)) / n)

  }

  # Result
  return(bw_select(loss = function(h) sapply(h, mise_star), optim = optim,
                   h_grid = h_grid, optim_par = optim_par,
                   optim_lower = optim_lower, optim_upper = optim_upper,
                   plot_it = plot_it, ylab = "MISE", name = "MISE"))

}

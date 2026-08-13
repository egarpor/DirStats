

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
#' # Sample
#' q <- 2
#' n <- 300
#' set.seed(42)
#' samp <- rbind(rotasym::r_vMF(n = n / 3, mu = c(rep(0, q), 1), kappa = 5),
#'               rotasym::r_vMF(n = n / 3, mu = c(rep(0, q), -1), kappa = 5),
#'               rotasym::r_vMF(n = n / 3, mu = c(1, rep(0, q)), kappa = 5))
#'
#' # Mixture fit
#' bic_vmf_mix(data = samp, plot_it = TRUE, verbose = TRUE)
#'
#' # Custom criterion: the AIC with the finite-sample correction
#' # 2 * k * (k + 1) / (n - k - 1), for k parameters and n observations
#' AICc <- function(object) {
#'   ll <- logLik(object)
#'   k <- attr(ll, "df")
#'   n <- attr(ll, "nobs")
#'   AIC(object) + ifelse(n - k - 1 > 0, 2 * k * (k + 1) / (n - k - 1), Inf)
#' }
#' bic_vmf_mix(data = samp, crit = AICc, plot_it = TRUE)$BICs
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
  get_BIC <- function(x) {
    if (is.na(x[1])) NA else crit(x)
  }

  # Screening of degenerate mixtures
  over_kappa_max <- function(x) {
    !is.na(x[1]) && any(norm2(x$theta) > kappa_max)
  }

  # Mixture fitting for 1:M_bound
  if (verbose) {

    cat("Fitting mixtures...\n")
    pb <- txtProgressBar(style = 3)

  }
  for (k in 1:M_bound) {

    # Fit mixture
    fit_mixs[[k]] <- tryCatch(movMF::movMF(x = data, k = k,
                                            maxiter = 100, start = "S"),
                              silent = TRUE, error = function(e) NA)
    BICs[k] <- get_BIC(fit_mixs[[k]])

    if (is.na(fit_mixs[[k]][1]) || over_kappa_max(fit_mixs[[k]])) {

      BICs[k] <- NA

    }

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

        # New estimation
        fit_mixs_new <- tryCatch(movMF::movMF(x = data, k = new_M_bound,
                                               maxiter = 100, start = "S"),
                                 silent = TRUE, error = function(e) NA)

        # Add to the previous estimation
        fit_mixs <- c(fit_mixs, list(fit_mixs_new))

        # Compute the new BICs
        BICs_new <- get_BIC(fit_mixs_new)
        if (is.na(fit_mixs_new[1]) || over_kappa_max(fit_mixs_new)) {

          BICs_new <- NA

        }
        BICs <- c(BICs, BICs_new)

        # Seek for the minimum. Note that NA and NaN values are discarded
        M_opt <- which.min(BICs)

        # Stop iteration if M_opt < new_M_bound - M_neig + 1; otherwise
        # increase the bound
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

    # Compute a more precise solution for the optimum
    fit_mixs_better <- tryCatch(movMF::movMF(x = data, k = M_opt,
                                             maxiter = 300, start = "S"),
                                silent = TRUE, error = function(e) NULL)

    # Take fit_mixs_better only if there are no convergence issues, it does not
    # worsen crit, and it is not degenerate, as a longer run may land on a
    # spiked mixture
      fit_mix <- fit_mixs[[M_opt]]
if (!is.null(fit_mixs_better) && !over_kappa_max(fit_mixs_better)) {

    BIC_better <- get_BIC(fit_mixs_better)
      if (!is.na(BIC_better) && BIC_better <= BICs[M_opt]) {

      fit_mix <- fit_mixs_better

}

    }

    # Extract parameters
    kappa_hat <- norm2(fit_mix$theta)
    mu_hat <- fit_mix$theta / kappa_hat
    p_hat <- fit_mix$alpha

  }

  # Display a graph?
  if (plot_it) {

    plot(seq_len(ifelse(iterative, new_M_bound, M_bound)), BICs, xlab = "M",
         ylab = crit_name, main = paste("Min", crit_name), pch = 19, cex = 1.5)
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
#' @param fit_mix output from \code{\link{bic_vmf_mix}}. Computed internally if
#' \code{NULL} (default).
#' @param mu,kappa,p mixture parameters. \code{mu} is the mean matrix of size
#' \code{c(length(p), q + 1)}, \code{kappa} is vector of \code{length(p)}
#' concentration parameters, and \code{p} is the vector of mixture proportions.
#' @param optim_par,optim_lower,optim_upper parameters passed to \code{par},
#' \code{lower}, and \code{upper} in \code{\link[stats]{optim}} when using the
#' \code{"L-BFGS-B"} method. Default to \code{0.25},
#' \code{sqrt(.Machine$double.eps)}, and \code{10}.
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
#' set.seed(42)
#' samp <- rotasym::r_vMF(n = n, mu = c(1, rep(0, q)), kappa = 2)
#'
#' # Mixture fit
#' fit_mix <- bic_vmf_mix(data = samp, plot_it = TRUE)
#'
#' # ROT
#' bw_dir_rot(samp)
#'
#' # AMI
#' bw_dir_ami(samp)
#' bw_dir_ami(samp, fit_mix = fit_mix)
#' bw_dir_ami(samp, fit_mix = fit_mix, L = function(x) exp(-x))
#'
#' # EMI
#' bw_dir_emi(samp)
#' bw_dir_emi(samp, fit_mix = fit_mix, optim = FALSE, plot_it = TRUE)
#' @name bw_dir_pi


#' @rdname bw_dir_pi
#' @export
bw_dir_rot <- function(data) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Set parameters
  n <- nrow(data)
  q <- ncol(data) - 1
  kappa <- kappa_ml(data = data)

  # Degenerate concentrations (uniform data or data collapsed at a single point)
  if (!is.finite(kappa) || kappa <= 0) {

    warning(paste0("Degenerate concentration estimate (kappa = ", kappa,
                   "); returning NA."))
    return(NA)

  }

  # Rule of thumb of Proposition 2 in Garcia-Portugues (2013), evaluated in log
  # scale for every q at once. Writing L(nu, x) = log(exp(-x) * I_nu(x)), the
  # exp(2 * kappa) factors of the numerator and denominator cancel, and the sum
  # in the denominator is accumulated with a log-sum-exp. Relying on
  # log_besselI_scaled() keeps the selector finite for arbitrarily large
  # concentrations. There is a typo in Garcia-Portugues (2013) in equation (6)
  # for q = 2, where it displays an incorrect extra kappa:
    # * WRONG:    kappa * (1 + 4 * kappa^2) * sinh(2 * kappa)
    # * CORRECT:          (1 + 4 * kappa^2) * sinh(2 * kappa)
# This is not relevant for the general case q >= 1 compted here
log_num <- log(4 * sqrt(pi)) +
    2 * log_besselI_scaled(nu = (q - 1) / 2, x = kappa)
  log_1 <- log(2 * q) +
    log_besselI_scaled(nu = (q + 1) / 2, x = 2 * kappa)
  log_2 <- log((2 + q) * kappa) +
    log_besselI_scaled(nu = (q + 3) / 2, x = 2 * kappa)
  log_max <- pmax(log_1, log_2)
  log_den <- ((q + 1) / 2) * log(kappa) + log(n) +
    log_max + log(exp(log_1 - log_max) + exp(log_2 - log_max))
  return(exp((log_num - log_den) / (q + 4)))

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
  if (length(p) != nrow(mu) || length(p) != length(kappa) ||
      length(kappa) != nrow(mu)) {

    stop("Check the dimension of the arguments")

  }

  # Psi operator for a single vMF(mu, kappa)
  Psi_vmf <- function(x, mu, kappa) {

    q <- length(mu) - 1
    log_const <- log_c_vmf(q = q, kappa = kappa)
    Psi <- drop(kappa * exp(log_const + kappa * x %*% mu) *
                  (-x %*% mu + kappa / q * (1 - (x %*% mu)^2)))
    return(Psi)

  }

  # Exact expression
  if (length(p) == 1) {

    # Evaluated in log scale (with a log-sum-exp for the two-term sum), so that
    # the curvature stays finite for large concentrations. The exp(2 * kappa)
    # factors of const and R cancel.
    log_const <- ((q + 1) / 2) * log(kappa) -
      ((q + 2) * log(2) + ((q + 1) / 2) * log(pi) + log(q) +
         2 * log_besselI_scaled(nu = (q - 1) / 2, x = kappa))
    log_1 <- log(2 * q) +
      log_besselI_scaled(nu = (q + 1) / 2, x = 2 * kappa)
    log_2 <- log((2 + q) * kappa) +
      log_besselI_scaled(nu = (q + 3) / 2, x = 2 * kappa)
    log_max <- pmax(log_1, log_2)
    R <- exp(log_const + log_max +
               log(exp(log_1 - log_max) + exp(log_2 - log_max)))

  # Numeric approximation
  } else {

    # Curvature integrand
    integrand <- function(x) {

      rowSums(sapply(seq_along(p), function(i)
        Psi_vmf(x = x, mu = mu[i, ], kappa = kappa[i])
      ) %*% p)^2

    }

    # Integration
    if (q == 1) {

      R <- tryCatch(int_cir(f = integrand, N = 500),
                    silent = TRUE, error = function(e) NULL)

    } else if (q == 2) {

      R <- tryCatch(int_sph(f = integrand),
                    silent = TRUE, error = function(e) NULL)

    } else if (q > 2) {

      R <- tryCatch(int_hypsph(f = integrand, q = q, M = 1e4),
                    silent = TRUE, error = function(e) NULL)

    }

    # To avoid numeric(0) result
    R <- ifelse(length(R) == 0, NA, R)

  }

  # Curvature term
  return(R)

}


#' @rdname bw_dir_pi
#' @export
bw_dir_emi <- function(data, fit_mix = NULL, optim = TRUE,
                       h_grid = exp(seq(log(0.05), log(1.5), l = 100)),
                       plot_it = TRUE, optim_par = 0.25,
                       optim_lower = sqrt(.Machine$double.eps),
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

  # Integrand for MISE function
  integrand <- function(x, h) {

    # Startup
    if (!is.matrix(x)) x <- t(x)
    q <- ncol(x) - 1

    # Expectation of the estimator, computed in log scale so that neither the
    # to avoid overflow/underflow for small h.
    h2 <- h^2
    log_inv_Cq_x <- sapply(seq_along(p), function(j) {
      nrm <- norm2(sweep(x / h2, 2, mu[j, ] * kappa[j], "+"))
      log(p[j]) + log_c_vmf(q = q, kappa = kappa[j]) -
        log_c_vmf(q = q, kappa = nrm)
    })
    if (is.null(dim(log_inv_Cq_x))) {

      log_inv_Cq_x <- matrix(log_inv_Cq_x, nrow = nrow(x), ncol = length(p))

    }

    # log(Eq) with the exponentially scaled Bessel
    log_Eq <- ((q + 1) / 2) * log(2 * pi) + (q - 1) * log(h) +
      log_besselI_scaled(nu = (q - 1) / 2, x = 1 / h2) + 1 / h2

    # Expectation of the estimator
    Ef <- rowSums(exp(log_inv_Cq_x - log_Eq))

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

  # Definition of mise_star
  mise_star <- function(h) {

    # Integral
    if (q == 1) {

      int <- int_cir(f = integrand, N = 500, h = h)

    } else if (q == 2) {

      int <- int_sph(f = integrand, h = h)

    } else {

      int <- int_hypsph(f = integrand, q = q, M = 1e4, h = h)

    }

    # Save time in computing quotients of ch. Done in log scale, and with the
    # exponentially scaled Bessel: the exp(1 / h^2) factors of Eq and of the
    # Bessel at 1 / h^2 cancel with the exp(2 / h^2) of the Bessel at 2 / h^2.
    log_b1 <- log_besselI_scaled(nu = (q - 1) / 2, x = 1 / h^2)
    log_b2 <- log_besselI_scaled(nu = (q - 1) / 2, x = 2 / h^2)
    log_Eq <- ((q + 1) / 2) * log(2 * pi) + (q - 1) * log(h) + log_b1
    log_Dq <- log_Eq + log_b1 - log_b2 + ((q - 1) / 2) * log(2)
    int <- int + exp(-log_Dq) / n
    return(int)

  }

  # Optimize versus grid search
  if (optim) {

    # Search for the minimum
    opt <- optim(fn = function(x) mise_star(h = x), par = optim_par,
                 lower = optim_lower, upper = optim_upper, method = "L-BFGS-B")

    # Result
    return(list(h_opt = opt$par, h_grid = NULL, MISE_opt = opt$value,
                MISE_grid = NULL))

  } else {

    # Grid search
    MISE <- sapply(h_grid, mise_star)

    # Index of the minimun CV
    ind_min <- which.min(MISE)

    # Warning if extreme
    if (h_grid[ind_min] == max(h_grid) || h_grid[ind_min] == min(h_grid)) {

      warning("h_opt at the extreme of h_grid")

    }

    # Plot?
    if (plot_it) {

      plot(h_grid, MISE, xlab = "h", ylab = "MISE", type = "o")
      rug(h_grid)
      abline(v = h_grid[ind_min], col = 2)

    }

    # Result
    return(list(h_opt = h_grid[ind_min], h_grid = h_grid,
                MISE_opt = MISE[ind_min], MISE_grid = MISE))

  }

}

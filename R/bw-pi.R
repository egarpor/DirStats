

#' @title Fitting mixtures of von Mises--Fisher distributions
#'
#' @description Fitting mixtures of von Mises--Fisher distributions by the
#' Expectation-Maximization algorithm, with determination of the optimal
#' number of mixture components.
#'
#' @inheritParams kde_dir
#' @param M_bound bound for the number of components in the mixtures. If it is
#' not enough, the search for the mixture with minimum \code{crit} will
#' continue from \code{M_bound + 1} if \code{iterative = TRUE}. Defaults to
#' \code{ceiling(log(nrow(data)))}.
#' @param M_neig number of neighbors explored around the optimal number
#' of mixture components. Defaults to \code{3}.
#' @param crit information criterion employed, either \code{"BIC"} (default),
#' \code{"AICc"} or \code{"AIC"}.
#' @param iterative keep exploring higher number of components if the optimum
#' is attained at \code{M_bound}? Defaults to \code{TRUE}.
#' @inheritParams bw_dir_cv
#' @param verbose display fitting progress? Defaults to \code{FALSE}.
#' @param kappa_max maximum value of allowed concentrations, to avoid numerical
#' instabilities. Defaults to \code{250}.
#' @return A list with entries:
#' \itemize{
#'   \item{\code{best_fit}: a list with estimated mixture parameters
#'   \code{mu_hat}, \code{kappa_hat}, and \code{p_hat} of the best-fitting
#'   mixture according to \code{crit}.}
#'   \item{\code{fit_mixs}: a list with of the fitted mixtures.}
#'   \item{\code{BICs}: a vector with the BICs (or other information criterion)
#'   of the fitted mixtures.}
#' }
#' @details
#' See Algorithm 3 in García-Portugués (2013). The Expectation-Maximization
#' fit is performed with \code{\link[movMF]{movMF}}.
#' @references
#' García-Portugués, E. (2013). Exact risk improvement of bandwidth selectors
#' for kernel density estimation with directional data. \emph{Electronic
#' Journal of Statistics}, 7:1655--1685.
#' \doi{10.1214/13-ejs821}
#'
#' Hornik, K. and Grün, B. (2014). movMF: An R Package for Fitting Mixtures of
#' von Mises--Fisher Distributions. \emph{Journal of Statistical Software},
#' 58(10):1--31. \doi{10.18637/jss.v058.i10}
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
#' @export
bic_vmf_mix <- function(data, M_bound = ceiling(log(nrow(data))), M_neig = 3,
                        crit = "BIC", iterative = TRUE, plot_it = FALSE,
                        verbose = FALSE, kappa_max = 250) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Create list for mixtures
  fit_mixs <- vector("list", length = M_bound)

  # Create vector for BICs
  BICs <- rep(NA, M_bound)

  # For computing the BICs
  get_BIC <- function(x) ifelse(is.na(x[1]), NA, get(crit)(x))

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

    # Avoid fitted mixtures with kappas > 250, they only cause problems
    # exp(250) = 3.746455e+108
    if (is.na(fit_mixs[[k]][1])) {

      BICs[k] <- NA

    } else if (any(norm2(fit_mixs[[k]]$theta) > kappa_max)) {

      BICs[k] <- NA

    }

    if (verbose) setTxtProgressBar(pb, k / M_bound)

  }

  # Seek for the minimum. Note that NA and NaN values are discarded
  M_opt <- which.min(BICs)

  # What if M_opt is near the maximum of the grid? Run the iterative procedure
  if (iterative) {

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
        if (is.na(fit_mixs_new[1])) {

          BICs_new <- NA

        } else if (any(norm2(fit_mixs_new$theta) > kappa_max)) {

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
    mu_hat <- mu_ml(data = data)
    kappa_hat <- kappa_ml(data = data)
    p_hat <- 1

  } else {

    # Compute a more precise solution for the optimum
    fit_mixs_better <- tryCatch(movMF::movMF(x = data, k = M_opt,
                                             maxiter = 300, start = "S"),
                                silent = TRUE, error = function(e) NULL)

    # Return fit_mixs_better if there are no convergence issues
    if (is.null(fit_mixs_better)) {

      fit_mix <- fit_mixs[[M_opt]]

    } else {

      fit_mix <- fit_mixs_better

    }

    # Extract parameters
    kappa_hat <- norm2(fit_mix$theta)
    mu_hat <- fit_mix$theta / kappa_hat
    p_hat <- fit_mix$alpha

  }

  # Display a graph?
  if (plot_it) {

    plot(1:ifelse(iterative, new_M_bound, M_bound), BICs, xlab = "M",
         ylab = crit, main = paste("Min", crit), pch = 19, cex = 1.5)
    points(M_opt, BICs[M_opt], pch = 19, cex = 1.5, col = 2)

  }

  # Final object
  return(list(best_fit =
                list(mu_hat = mu_hat, kappa_hat = kappa_hat, p_hat = p_hat),
              fits = fit_mixs, BICs = BICs))

}


#' @title Plug-in bandwidth selectors for directional data
#'
#' @description Plug-in bandwidth selectors for kernel density estimation
#' with directional data, including Rule-Of-Thumb (ROT),
#' Asymptotic MIxtures (AMI), and Exact MIxtures (EMI).
#'
#' @inheritParams kde_dir
#' @inheritParams bic_vmf_mix
#' @inheritParams bw_dir_cv
#' @param fit_mix output from \code{\link{bic_vmf_mix}}. Computed internally
#' if \code{NULL} (default).
#' @param mu,kappa,p mixture parameters. \code{mu} is the mean matrix of size
#' \code{c(length(p), q + 1)}, \code{kappa} is vector of \code{length(p)}
#' concentration parameters, and \code{p} is the vector of mixture proportions.
#' @param optim_par,optim_lower,optim_upper parameters passed to \code{par},
#' \code{lower}, and \code{upper} in \code{\link[stats]{optim}} when using
#' the \code{"L-BFGS-B"} method. Default to \code{0.25}, \code{0.06}
#' (to avoid numerical instabilities), and \code{10}.
#' @return Selected bandwidth for \code{bw_dir_rot} and \code{bw_dir_ami}.
#' \code{bw_dir_emi} returns a list with entries:
#' \itemize{
#'   \item{\code{h_opt}: selected bandwidth.}
#'   \item{\code{h_grid}: \code{h_grid}, if used (otherwise \code{NULL}).}
#'   \item{\code{MISE_opt}: minimum of the MISE loss.}
#'   \item{\code{MISE_grid}: value of the MISE function at \code{h_grid}, if
#'    used (otherwise \code{NULL}).}
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
#' \doi{10.1214/13-ejs821}
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

  # Cases for different dimensions
  if (q == 1) {

    num <- 4 * sqrt(pi) * besselI(nu = 0, x = kappa)^2
    den <- kappa * (2 * besselI(nu = 1, x = 2 * kappa) +
                      3 * kappa * besselI(nu = 2, x = 2 * kappa)) * n

  } else if (q == 2) {

    e <- exp(-2 * kappa)
    e2 <- e^2
    sinh_scaled <- 0.5 * (1 - e)
    cosh_scaled2 <- 0.5 * (1 + e2)
    sinh_scaled2 <- 0.5 * (1 - e2)
    num <- 8 * sinh_scaled^2
    den <- kappa * (-2 * kappa * cosh_scaled2 +
                      (1 + 4 * kappa^2) * sinh_scaled2) * n
    # Caution! Typo in G-P (2013) in equation (6) for q = 2. It displays an
    # incorrect extra kappa:
    # * WRONG:    kappa * (1 + 4 * kappa^2) * sinh(2 * kappa)
    # * CORRECT:          (1 + 4 * kappa^2) * sinh(2 * kappa)

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
  if (length(p) != nrow(mu) | length(p) != length(kappa) |
      length(kappa) != nrow(mu)) {

    stop("Check the dimension of the arguments")

  }

  # Psi operator for a single vMF(mu, kappa)
  Psi_vmf <- function(x, mu, kappa) {

    q <- length(mu) - 1
    log_const <- rotasym::c_vMF(p = q + 1, kappa = kappa, log = TRUE)
    Psi <- drop(kappa * exp(log_const + kappa * x %*% mu) *
                  (-x %*% mu + kappa / q * (1 - (x %*% mu)^2)))
    return(Psi)

  }

  # Exact expression
  if (length(p) == 1) {

    const <- kappa^((q + 1) / 2) /
      (2^(q + 2) * pi^((q + 1) / 2) *
         besselI(nu = (q - 1) / 2, x = kappa, expon.scaled = TRUE)^2 * q)
    R <- 2 * q * besselI(nu = (q + 1) / 2, x = 2 * kappa,
                         expon.scaled = TRUE) +
      (2 + q) * kappa * besselI(nu = (q + 3) / 2, x = 2 * kappa,
                                expon.scaled = TRUE)
    R <- const * R

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
                       plot_it = TRUE, optim_par = 0.25, optim_lower = 0.06,
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

    # Constants depending on x and h
    h2 <- h^2
    inv_Cq_x <- sapply(seq_along(p), function(j) {
      p[j] * rotasym::c_vMF(p = q + 1, kappa = kappa[j]) /
        rotasym::c_vMF(p = q + 1, kappa =
                         norm2(sweep(x / h2, 2, mu[j, ] * kappa[j], "+")))
    })
    if (is.null(dim(inv_Cq_x))) {

      inv_Cq_x <- matrix(inv_Cq_x, nrow = nrow(x), ncol = length(p))

    }

    # Save time in computing quotients of ch
    bb <- besselI(x = 1 / h2, nu = (q - 1) / 2)
    Eq <- (2 * pi)^((q + 1) / 2) * h^(q - 1) * bb

    # Expectation of the estimator
    Ef <- rowSums(inv_Cq_x) / Eq

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

    # Save time in computing quotients of ch
    bb <- besselI(x = 1 / h^2, nu = (q - 1) / 2, expon.scaled = TRUE)
    Eq <- (2 * pi)^((q + 1) / 2) * h^(q - 1) * bb
    Dq <- Eq * bb /
      besselI(x = 2 / h^2, nu = (q - 1) / 2, expon.scaled = TRUE) *
      2^((q - 1) / 2)
    int <- int + 1 / (Dq * n)
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
    if (h_grid[ind_min] == max(h_grid) | h_grid[ind_min] == min(h_grid)) {

      warning("h_opt at the extreme of h_grid")

    }

    # Plot?
    if (plot_it) {

      plot(h_grid, MISE, type = "l", xlab = "h", ylab = "MISE")
      rug(h_grid)
      abline(v = h_grid[ind_min], col = 2)

    }

    # Result
    return(list(h_opt = h_grid[ind_min], h_grid = h_grid,
                MISE_opt = MISE[ind_min], MISE_grid = MISE))

  }

}

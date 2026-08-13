

#' @title Stable evaluation of scaled log-Bessel functions and the von
#' Mises--Fisher normalizing constant
#'
#' @description Numerically stable evaluation of the logarithm of the
#' exponentially scaled modified Bessel function of the first kind,
#' \eqn{\log(e^{-x} \mathcal{I}_\nu(x))}, and of the (log) von Mises--Fisher
#' normalizing constant that builds on it. Both remain finite and accurate for
#' large orders \eqn{\nu} and large arguments \eqn{x}, where a direct call to
#' \code{\link[base]{besselI}} would overflow, underflow, or lose precision.
#'
#' @param nu order of the Bessel function. A scalar or a vector conformable
#' with \code{x}.
#' @param x vector of non-negative evaluation points.
#' @param p dimension of the ambient space \eqn{R^p} (so the data lies on
#' \eqn{S^{p-1}}).
#' @param kappa vector of non-negative concentration parameters.
#' @param spline if \code{TRUE} (default), evaluate the scaled log-Bessel with a
#' fast cubic-spline interpolation for orders \code{nu} in
#' \code{seq(0, 50, by = 0.5)}, and fall back to \code{\link[base]{besselI}} for
#' any other order. If \code{FALSE}, always use \code{\link[base]{besselI}}
#' together with the asymptotic expansions. The spline is only effective when
#' \code{nu} is a single value on that grid.
#' @details \code{log_besselI_scaled} combines four regimes so that the result
#' is finite and accurate over the whole \eqn{(\nu, x)} plane:
#' \itemize{
#'   \item \eqn{x < 10^{-5}}: the exact leading small-argument series
#'   \eqn{\nu \log(x / 2) - \log\Gamma(\nu + 1) - x}.
#'   \item \eqn{x \ge 10^{4}}: the large-argument asymptotic expansion
#'   \code{\link[Bessel]{besselIasym}}.
#'   \item \eqn{\nu > 50}: the large-order (DLMF 10.41) asymptotic expansion
#'   \code{\link[Bessel]{besselI.nuAsym}}, which meets the
#'   \eqn{\nu \le 50} regime smoothly.
#'   \item otherwise: \code{\link[base]{besselI}} with
#'   \code{expon.scaled = TRUE}, or, when \code{spline = TRUE}, the interpolated
#'   spline. Any underflow of \code{besselI} at small \eqn{x} and large
#'   \eqn{\nu} is filled with the small-argument series.
#' }
#' The cutoff \eqn{\nu = 50} is where \code{besselI.nuAsym(k.max = 5)} is
#' accurate to below \eqn{10^{-10}}, so it takes over from the tabulated orders
#' with a smooth join.
#'
#' \code{fast_log_c_vMF} returns
#' \eqn{\log C_p(\kappa)}, the logarithm of the von Mises--Fisher normalizing
#' constant \eqn{C_p(\kappa) = \kappa^{p / 2 - 1} / ((2\pi)^{p / 2}
#' \mathcal{I}_{p / 2 - 1}(\kappa))}, evaluated stably through
#' \code{log_besselI_scaled}.
#' @return
#' \itemize{
#'   \item \code{log_besselI_scaled}: the evaluation of
#'   \eqn{\log(e^{-x} \mathcal{I}_\nu(x))}, a vector conformable with \code{nu}
#'   and \code{x}.
#'   \item \code{fast_log_c_vMF}: the evaluation of \eqn{\log C_p(\kappa)}, a
#'   vector of size \code{length(kappa)}.
#' }
#' @references
#' DLMF (2024). \emph{NIST Digital Library of Mathematical Functions}. Section
#' 10.41. \url{https://dlmf.nist.gov/10.41}
#' @examples
#' # Scaled log-Bessel across regimes
#' log_besselI_scaled(nu = 0.5, x = c(1e-6, 1, 1e5))
#' # The spline (default) and the direct evaluation agree closely
#' log_besselI_scaled(nu = 0.5, x = 2)
#' log_besselI_scaled(nu = 0.5, x = 2, spline = FALSE)
#'
#' # Stable where besselI() is not
#' log_besselI_scaled(nu = 400, x = 20)
#'
#' # (Log) von Mises--Fisher normalizing constant, uniform case included
#' fast_log_c_vMF(p = 3, kappa = c(0, 1, 10))
#' @name bessel
NULL


#' @rdname bessel
#' @export
log_besselI_scaled <- function(nu, x, spline = TRUE) {

  # Conformability of nu and x
  if (!(length(x) == 1 || length(nu) == 1 || length(nu) == length(x))) {

    stop("nu and x must have equal lengths, or one of them must be a scalar.")

  }
  if (anyNA(nu)) {

    stop("nu must not contain NAs.")

  }

  # Zero-length input gives zero-length output, as in besselI()
  if (length(x) == 0 || length(nu) == 0) {

    return(numeric(0))

  }

  # Broadcast to a common length; NA positions in x propagate to the output
  n <- max(length(x), length(nu))
  x <- rep_len(x, n)
  nu <- rep_len(nu, n)
  out <- rep(NA_real_, n)
  ok <- !is.na(x)
  if (!any(ok)) {

    return(out)

  }
  xo <- x[ok]
  nuo <- nu[ok]
  res <- rep(NA_real_, length(xo))

  # Largest tabulated order; beyond it the large-order asymptotic takes over,
  # so spline = TRUE meets the asymptotic regime smoothly
  nu_max <- max(bessel_nus)

  # Spline path: requires a single tabulated order (the von Mises--Fisher orders
  # (p - 2) / 2 are half-integers, hence on the grid); otherwise fall back to
  # the direct/asymptotic evaluation
  if (spline && all(nuo == nuo[1L]) && nuo[1L] <= nu_max &&
        any(abs(bessel_nus - nuo[1L]) < 1e-8)) {

    nu1 <- nuo[1L]
    ind_small <- xo < bessel_x_lo
    ind_hi <- xo >= bessel_x_hi
    ind_spl <- !ind_small & !ind_hi

    # Cubic spline in u = log(x) over [x_lo, x_hi)
    if (any(ind_spl)) {

      key <- sprintf("%03d", round(10 * nu1))
      res[ind_spl] <- get(paste0("log_besselI_scaled_spline_", key))(
        log(xo[ind_spl]))

    }

    # Small-argument series below x_lo (exact leading term; the nu = 0 branch
    # avoids 0 * log(0) = NaN)
    if (any(ind_small)) {

      xs <- xo[ind_small]
      res[ind_small] <- if (nu1 == 0) -xs else
        nu1 * log(xs / 2) - lgamma(nu1 + 1) - xs

    }

    # Large-argument asymptotic at/above x_hi
    if (any(ind_hi)) {

      res[ind_hi] <- Bessel::besselIasym(x = xo[ind_hi], nu = nu1, k.max = 10,
                                         expon.scaled = TRUE, log = TRUE)

    }

  } else {

    # Direct evaluation with asymptotic fallbacks. The large-order expansion
    # (DLMF 10.41) takes precedence over the large-argument one (the latter is
    # only valid for x >> nu^2)
    ind_nu <- nuo > nu_max
    ind_x <- !ind_nu & (xo >= bessel_x_hi)
    ind_std <- !ind_nu & !ind_x

    if (any(ind_nu)) {

      res[ind_nu] <- Bessel::besselI.nuAsym(x = xo[ind_nu], nu = nuo[ind_nu],
                                            k.max = 5, expon.scaled = TRUE,
                                            log = TRUE)

    }
    if (any(ind_x)) {

      res[ind_x] <- Bessel::besselIasym(x = xo[ind_x], nu = nuo[ind_x],
                                        k.max = 10, expon.scaled = TRUE,
                                        log = TRUE)

    }
    if (any(ind_std)) {

      v <- log(besselI(x = xo[ind_std], nu = nuo[ind_std],
                       expon.scaled = TRUE))

      # Fill any underflow (large order + small x) with the exact leading
      # small-argument series, so the result is finite for all tabulated orders
      bad <- !is.finite(v)
      if (any(bad)) {

        xs <- xo[ind_std][bad]
        ns <- nuo[ind_std][bad]
        v[bad] <- ifelse(ns == 0, -xs, ns * log(xs / 2) - lgamma(ns + 1) - xs)

      }
      res[ind_std] <- v

    }

  }

  out[ok] <- res
  return(out)

}


#' @rdname bessel
#' @export
fast_log_c_vMF <- function(p, kappa, spline = TRUE) {

  log_bessel <- log_besselI_scaled(nu = 0.5 * (p - 2), x = kappa,
                                   spline = spline)
  log_c_vMF <- (0.5 * (p - 2)) * log(kappa) - (0.5 * p) * log(2 * pi) -
    kappa - log_bessel

  # Uniform density limit (kappa = 0): the reciprocal surface area of S^{p - 1}
  log_c_vMF[kappa == 0] <- -rotasym::w_p(p = p, log = TRUE)
  return(log_c_vMF)

}

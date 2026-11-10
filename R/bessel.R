

#' @title Stable evaluation of the scaled log-Bessel function
#'
#' @description Numerically stable evaluation of the logarithm of the
#' exponentially scaled modified Bessel function of the first kind,
#' \eqn{\log(e^{-x} \mathcal{I}_\nu(x))}. It remains finite and accurate for
#' large orders \eqn{\nu} and large arguments \eqn{x}, where a direct call to
#' \code{\link[base]{besselI}} would overflow, underflow, or lose precision.
#'
#' @param nu order of the Bessel function. A scalar or a vector conformable
#' with \code{x}.
#' @param x vector of non-negative evaluation points.
#' @param spline if \code{TRUE} (default), evaluate the scaled log-Bessel with a
#' fast cubic-spline interpolation for orders \code{nu} in
#' \code{seq(0, 50, by = 0.5)}, and fall back to \code{\link[base]{besselI}} for
#' any other order. If \code{FALSE}, always use \code{\link[base]{besselI}}
#' together with the asymptotic expansions. The spline is only effective when
#' \code{nu} is a single value on that grid.
#' @return The evaluation of \eqn{\log(e^{-x} \mathcal{I}_\nu(x))}, a vector
#' conformable with \code{nu} and \code{x}.
#' @details The function combines four regimes so that the result is finite and
#' accurate over the whole \eqn{(\nu, x)} plane:
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
#' @references
#' DLMF (2024). \emph{NIST Digital Library of Mathematical Functions}. Section
#' 10.41. \url{https://dlmf.nist.gov/10.41}
#' @seealso \code{\link{log_c_vmf}}, \code{\link[base]{besselI}}.
#' @examples
#' # Scaled log-Bessel across regimes
#' log_besselI_scaled(nu = 0.5, x = c(1e-6, 1, 1e5))
#' # The spline (default) and the direct evaluation agree closely
#' log_besselI_scaled(nu = 0.5, x = 2)
#' log_besselI_scaled(nu = 0.5, x = 2, spline = FALSE)
#'
#' # Stable where besselI() is not
#' log_besselI_scaled(nu = 400, x = 20)
#' @export
log_besselI_scaled <- function(nu, x, spline = TRUE) {

  # Conformability of nu and x
  if (!(length(x) == 1 || length(nu) == 1 || length(nu) == length(x))) {

    stop("nu and x must have equal lengths, or one of them must be a scalar.")

  }
  if (anyNA(nu)) {

    stop("nu must not contain NAs.")

  }

  # Largest tabulated order, beyond which the large-order asymptotic takes over
  nu_max <- max(bessel_nus)

  # Fast path for a single tabulated order with no NAs and all x in
  # [bessel_x_lo, bessel_x_hi) skipping the partitioning below.
  if (spline && length(nu) == 1 && length(x) > 0 && !anyNA(x) &&
      min(x) >= bessel_x_lo && max(x) < bessel_x_hi && nu <= nu_max &&
      any(abs(bessel_nus - nu) < 1e-8)) {

    return(bessel_splines[[round(2 * nu) + 1]](log(x)))

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

  # Spline path: requires a single tabulated order (the von Mises--Fisher orders
  # (p - 2) / 2 are half-integers, hence on the grid); otherwise fall back to
  # the direct/asymptotic evaluation.
  if (spline && all(nuo == nuo[1L]) && nuo[1L] <= nu_max &&
        any(abs(bessel_nus - nuo[1L]) < 1e-8)) {

    nu1 <- nuo[1L]
    ind_small <- xo < bessel_x_lo
    ind_hi <- xo >= bessel_x_hi
    ind_spl <- !ind_small & !ind_hi

    # Cubic spline in u = log(x) over [x_lo, x_hi)
    if (any(ind_spl)) {

      res[ind_spl] <- bessel_splines[[round(2 * nu1) + 1]](log(xo[ind_spl]))

    }

    # Small-argument series below x_lo
    if (any(ind_small)) {

      res[ind_small] <- log_besselI_series(x = xo[ind_small], nu = nu1)

    }

    # Large-argument asymptotic at/above x_hi
    if (any(ind_hi)) {

      res[ind_hi] <- Bessel::besselIasym(x = xo[ind_hi], nu = nu1, k.max = 10,
                                         expon.scaled = TRUE, log = TRUE)

    }

  } else {

    # Direct evaluation with asymptotic fallbacks. The large-order expansion
    # (DLMF 10.41) takes precedence over the large-argument one (the latter is
    # only valid for x >> nu^2).
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
      # small-argument series, so the result is finite for all tabulated orders.
      bad <- !is.finite(v)
      if (any(bad)) {

        v[bad] <- log_besselI_series(x = xo[ind_std][bad],
                                     nu = nuo[ind_std][bad])

      }
      res[ind_std] <- v

    }

  }

  out[ok] <- res
  return(out)

}


#' @title Small-argument series of the scaled log-Bessel function
#'
#' @description Leading term of \eqn{\log(e^{-x} \mathcal{I}_\nu(x))} as
#' \eqn{x \to 0}, \eqn{\nu \log(x / 2) - \log\Gamma(\nu + 1) - x}.
#'
#' @inheritParams log_besselI_scaled
#' @return A vector of length \code{length(x)}.
#' @noRd
log_besselI_series <- function(x, nu) {

  nu <- rep_len(nu, length(x))
  res <- nu * log(x / 2) - lgamma(nu + 1) - x
  ind_0 <- nu == 0
  res[ind_0] <- -x[ind_0]
  return(res)

}

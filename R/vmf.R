

#' @title Von Mises--Fisher distribution utilities
#'
#' @description Maximum likelihood estimation for the von Mises--Fisher
#' distribution and evaluation of density mixtures.
#'
#' @inheritParams kde_dir
#' @param mu,kappa,p mixture parameters. \code{mu} is the mean matrix of size
#' \code{c(length(p), q + 1)}, \code{kappa} is vector of \code{length(p)}
#' concentration parameters, and \code{p} is the vector of mixture proportions.
#' @param norm enforce normalization of \code{x} internally? Defaults to
#' \code{FALSE}.
#' @param min_kappa,max_kappa initial search interval (in \code{kappa}) for the
#' maximum likelihood estimate. The interval is automatically extended if the
#' root lies outside. Default to \code{1e-4} and \code{1e4}.
#' @param R resultant length, a scalar. If \code{NULL} (default), it is computed
#' from \code{data}.
#' @param ... further parameters passed to \code{\link{uniroot}}.
#' @return
#' \itemize{
#'   \item \code{mu_ml}: the estimated mean direction, a vector of size
#'   \code{q + 1}.
#'   \item \code{kappa_ml}: the estimated concentration parameter, a scalar.
#'   \item \code{d_mixvmf}: the evaluated mixture density, a vector of size
#'   \code{nx}.
#' }
#' @examples
#' # Sample
#' n <- 50
#' q <- 2
#' samp <- rotasym::r_vMF(n = n, mu = c(1, rep(0, q)), kappa = 2)
#'
#' # Estimates
#' mu_ml(samp)
#' kappa_ml(samp)
#'
#' # Mixture
#' x <- to_cir(seq(0, 2 * pi, l = 200))
#' dens <- d_mixvmf(x = x, mu = rbind(c(-1, 0), c(0, 1), c(1, 0)),
#'                  kappa = 1:3, p = c(0.5, 0.2, 0.3))
#' plot(to_rad(x), dens, type = "l")
#' @seealso \code{\link{kde_dir}}, which uses the von Mises--Fisher kernel, and
#' \code{\link{bic_vmf_mix}} for fitting mixtures of von Mises--Fisher
#' distributions.
#' @name vmf


#' @rdname vmf
#' @export
kappa_ml <- function(data, min_kappa = 1e-4, max_kappa = 1e4, R = NULL, ...) {

  # Compute the resultant length if not provided
  p <- ncol(data)
  if (is.null(R)) {

    R <- norm2(colMeans(data, na.rm = TRUE))

  }

  # Degenerate resultant lengths: R = 0 (uniform data) means no concentration,
  # R = 1 (all data identical) means an infinite one
  if (R <= 0) {

    return(0)

  }
  if (R >= 1) {

    return(Inf)

  }

  # Solve A_p(kappa) = I_{p / 2}(kappa) / I_{p / 2 - 1}(kappa) = R for kappa
  # with the root is sought over log(kappa) with an auto-extending search
  # interval.
  A_p <- function(log_kappa) {

    kappa <- exp(log_kappa)
    exp(log_besselI_scaled(nu = p / 2, x = kappa) -
          log_besselI_scaled(nu = p / 2 - 1, x = kappa)) - R

  }
  root <- uniroot(A_p, lower = log(min_kappa), upper = log(max_kappa),
                  extendInt = "upX", ...)$root
  return(exp(root))

}


#' @rdname vmf
#' @export
mu_ml <- function(data) {

  return(normalize(colSums(data, na.rm = TRUE)))

}


#' @rdname vmf
#' @export
d_mixvmf <- function(x, mu, kappa, p, norm = FALSE) {

  # Dimension checks
  if (!is.matrix(x)) x <- t(x)
  if (!is.matrix(mu)) mu <- t(mu)
  q <- ncol(x) - 1
  stopifnot(ncol(mu) == q + 1)

  # Mixture components checks
  if (length(p) != length(kappa) || length(kappa) != nrow(mu) ||
      length(p) != nrow(mu)) {

    stop("Check mixtures arguments")

  }

  # Normalize
  if (norm) x <- normalize(x)

  # Density components
  dens <- sapply(seq_along(p), function(i) {
    p[i] * exp(rotasym::c_vMF(p = q + 1, kappa = kappa[i], log = TRUE) +
                 kappa[i] * x %*% mu[i, ])
  })

  # Add components
  if (!is.matrix(dens)) dens <- t(dens)
  return(rowSums(dens))

}

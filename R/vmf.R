

#' @title Von Mises--Fisher distribution utilities
#'
#' @description Maximum likelihood estimation for the von Mises--Fisher
#' distribution and evaluation of density mixtures.
#'
#' @inheritParams kde_dir
#' @param mu,kappa,p mixture parameters. \code{mu} is the mean matrix of size
#' \code{c(length(p), q + 1)}, \code{kappa} is vector of \code{length(p)}
#' concentration parameters, and \code{p} is the vector of mixture proportions.
#' @param norm enforce normalization of \code{x} internally? Defaults
#' to \code{FALSE}.
#' @return Estimated vector mean (\code{mu_ml}) or concentration parameter
#' (\code{kappa_ml}). A vector of length \code{nx} for \code{d_mixvmf}.
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
#' @name vmf


#' @rdname vmf
#' @export
kappa_ml <- function(data) {

  p <- ncol(data)
  R <- norm2(colMeans(data, na.rm = TRUE))
  A_p <- function(kappa) besselI(x = kappa, nu = p / 2) /
    besselI(x = kappa, nu = p / 2 - 1) - R
  return(uniroot(A_p, lower = 1e-04, upper = 100, tol = 1e-4)$root)

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
  if (length(p) != length(kappa) | length(kappa) != nrow(mu) |
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


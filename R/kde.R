

#' @title Directional kernel density estimator
#'
#' @description Kernel density estimation with directional data as in
#' the estimator of Bai et al. (1988).
#'
#' @param x evaluation points, a matrix of size \code{c(nx, q + 1)}.
#' @param data directional data, a matrix of size \code{c(n, q + 1)}.
#' @param h bandwidth, a scalar for \code{kde_dir}. Can be a vector
#' for \code{c_h}.
#' @param L kernel function. Set internally to \code{function(x) exp(-x)}
#' (von Mises--Fisher kernel) if \code{NULL} (default).
#' @param q dimension of \eqn{S^q}, \eqn{q\ge 1}.
#' @return \code{kde_dir} returns a vector of size \code{nx} with the
#' evaluated kernel density estimator. \code{c_h} returns the normalizing
#' constant for the kernel, a vector of length \code{length(h)}.
#' \code{lambda_L}, \code{b_L}, and \code{d_L} return moments of \code{L}.
#' @details
#' \code{data} is not checked to have unit norm, so the user must be careful.
#' When \code{L = NULL}, faster FORTRAN code is employed.
#' @references
#' Bai, Z. D., Rao, C. R., and Zhao, L. C. (1988). Kernel estimators of
#' density function of directional data. \emph{Journal of Multivariate
#' Analysis}, 27(1):24--39.
#' \url{https://doi.org/10.1016/0047-259X(88)90113-3}
#' @examples
#' # Sample
#' n <- 50
#' q <- 3
#' samp <- rotasym::r_vMF(n = n, mu = c(1, rep(0, q)), kappa = 2)
#'
#' # Evaluation points
#' x <- rbind(diag(1, nrow = q + 1), diag(-1, nrow = q + 1))
#'
#' # kde_dir
#' kde_dir(x = x, data = samp, h = 0.5, L = NULL)
#' kde_dir(x = x, data = samp, h = 0.5, L = function(x) exp(-x))
#'
#' # c_h
#' c_h(h = 0.5, q = q, L = NULL)
#' c_h(h = 0.5, q = q, L = function(x) exp(-x))
#'
#' # b_L
#' b_L(L = NULL, q = q)
#' b_L(L = function(x) exp(-x), q = q)
#'
#' # d_L
#' d_L(L = NULL, q = q)
#' d_L(L = function(x) exp(-x), q = q)
#'
#' # lambda_L
#' lambda_L(L = NULL, q = q)
#' lambda_L(L = function(x) exp(-x), q = q)
#' @useDynLib DirStats
#' @export
kde_dir <- function(x, data, h, L = NULL) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Dimensions of data
  n <- nrow(data)
  q <- ncol(data) - 1

  # Transform x into matrix if vector
  if (is.null(dim(x))) {

    x <- matrix(x, ncol = length(x), byrow = TRUE)

  }
  nx <- nrow(x)
  stopifnot(ncol(x) == q + 1)

  # von Mises--Fisher kernel?
  if (is.null(L)) {

    is_vm <- TRUE

  } else {

    is_vm <- FALSE
    ch <- c_h(h = h, q = q, L = L)

  }
  stopifnot(length(h) == 1)

  # Normalizing constant
  ch <- c_h(h = h, q = q, L = switch(is_vm + 1, L, NULL))

  # Approach type
  if (is_vm) {

    # Ensure data types for kde_dir
    x <- matrix(as.double(x), nrow = nx, ncol = q + 1)
    data_dir <- matrix(as.double(data), nrow = n, ncol = q + 1)

    # Call to kde_dir
    kde <- .Fortran("kde_dir_vmf", x = x, daa_dir = data_dir,
                    h = as.double(h), ch = as.double(ch),
                    n = as.integer(n), q = as.integer(q),
                    nx = as.integer(nx), kde = double(nx))$kde

  } else {

    kde <- (1 - data %*% t(x)) / h^2
    kde <- ch * colMeans(L(kde))

  }

  # Result
  return(kde)

}


#' @rdname kde_dir
#' @export
c_h <- function(h, q, L = NULL) {

  # von Mises--Fisher kernel?
  if (is.null(L)) {

    # Analytical expression
    ch <- (2 * pi)^((q + 1)/2) * h^(q - 1) *
      besselI(x = 1/h^2, nu = (q - 1)/2, expon.scaled = TRUE)

  } else {

    # Log-area of Omega_{q - 1}
    w_q <- exp(log(2) + q / 2 * log(pi) - lgamma(q / 2))

    # Numerical integration
    ch <- sapply(h, function(hh) {
      integrate(function(t) {
        L((1 - t) / hh^2) * (1 - t^2)^(q / 2 - 1)
        }, lower = -1 + 1e-5, upper = 1 - 1e-5, rel.tol = 1e-5,
        subdivisions = 200, stop.on.error = FALSE)$value * w_q
      })

  }

  # Normalizing constant
  return(1 / ch)

}


#' @rdname kde_dir
#' @export
lambda_L <- function(L = NULL, q) {

  # von Mises--Fisher kernel?
  if (is.null(L)) {

    # Analytical expression
    lambda <- (2 * pi)^(q / 2)

  } else {

    # Numerical integration
    lambda <- 2^(q / 2 - 1) * 2 * pi^(q / 2) / gamma(q / 2) *
      integrate(function(r) L(r) * r^(q / 2 - 1), lower = 0, upper = 100,
                rel.tol = 1e-5, subdivisions = 200, stop.on.error = FALSE)$value

  }

  # Moment
  return(lambda)

}


#' @rdname kde_dir
#' @export
b_L <- function(L = NULL, q) {

  # von Mises--Fisher kernel?
  if (is.null(L)) {

    # Analytical expression
    b <- q / 2

  } else {

    # Numerical integration
    b <- integrate(function(r) L(r) * r^(q / 2), lower = 0, upper = 100,
                   rel.tol = 1e-5, subdivisions = 200,
                   stop.on.error = FALSE)$value /
      integrate(function(r) L(r) * r^(q / 2 - 1), lower = 0, upper = 100,
                rel.tol = 1e-5, subdivisions = 200,
                stop.on.error = FALSE)$value

  }

  # Moment
  return(b)

}


#' @rdname kde_dir
#' @export
d_L <- function(L = NULL, q) {

  # von Mises--Fisher kernel?
  if (is.null(L)) {

    # Analytical expression
    d <- 2^(-q / 2)

  } else {

    # Numerical integration
    d <- integrate(function(r) L(r)^2 * r^(q / 2 - 1), lower = 0, upper = 100,
                   rel.tol = 1e-5, subdivisions = 200,
                   stop.on.error = FALSE)$value /
      integrate(function(r) L(r) * r^(q / 2 - 1), lower = 0, upper = 100,
                rel.tol = 1e-5, subdivisions = 200,
                stop.on.error = FALSE)$value

  }

  # Moment
  return(d)

}


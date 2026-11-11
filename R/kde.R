

#' @title Directional kernel density estimator
#'
#' @description Kernel density estimation with directional data on
#' \eqn{\mathbb{S}^q} as in the estimator of Bai et al. (1988).
#'
#' @param x evaluation points on \eqn{\mathbb{S}^q}, a matrix of size
#' \code{c(nx, q + 1)}.
#' @param data directional data on \eqn{\mathbb{S}^q}, a matrix of size
#' \code{c(n, q + 1)}.
#' @param h bandwidth, a scalar for \code{kde_dir}. Can be a vector for
#' \code{c_h}.
#' @param L kernel function. Set internally to \code{function(x) exp(-x)} (von
#' Mises--Fisher kernel) if \code{NULL} (default).
#' @param q dimension of \eqn{\mathbb{S}^q}, \eqn{q\ge 1}.
#' @param log return the normalizing constant in log scale? Defaults to
#' \code{FALSE}.
#' @return
#' \itemize{
#'   \item \code{kde_dir}: the evaluated kernel density estimator, a vector of
#'   size \code{nx}.
#'   \item \code{c_h}: the normalizing constant of the kernel, a vector of size
#'   \code{length(h)}.
#' }
#' @details
#' \code{data} is not checked to have unit norm, so the user must be careful.
#' @references
#' Bai, Z. D., Rao, C. R., and Zhao, L. C. (1988). Kernel estimators of
#' density function of directional data. \emph{Journal of Multivariate
#' Analysis}, 27(1):24--39.
#' \doi{10.1016/0047-259X(88)90113-3}
#' @seealso \code{\link{bw_dir_cv}} and \code{\link{bw_dir_pi}} for selecting
#' the bandwidth \code{h}, and \code{\link{vmf}} for the underlying von
#' Mises--Fisher utilities.
#' @examples
#' # Sample
#' n <- 50
#' q <- 3
#' if (requireNamespace("rotasym", quietly = TRUE)) {
#'
#'   set.seed(42)
#'   samp <- rotasym::r_vMF(n = n, mu = c(1, rep(0, q)), kappa = 2)
#'
#'   # Evaluation points
#'   x <- rbind(diag(1, nrow = q + 1), diag(-1, nrow = q + 1))
#'
#'   # kde_dir
#'   print(kde_dir(x = x, data = samp, h = 0.5, L = NULL))
#'   print(kde_dir(x = x, data = samp, h = 0.5, L = function(x) exp(-x)))
#'
#' }
#'
#' # c_h
#' c_h(h = 0.5, q = q, L = NULL)
#' c_h(h = 0.5, q = q, L = function(x) exp(-x))
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
  stopifnot(length(h) == 1)

  # Normalizing constant (L = NULL triggers the closed-form vMF expression)
  ch <- c_h(h = h, q = q, L = L)

  # Von Mises--Fisher kernel?
  if (is.null(L)) {

    L <- function(z) exp(-z)

  }

  # Kernel sums over blocks of evaluation points, so the memory stays flat in n
  kde <- numeric(nx)
  for (ind in row_blocks(n = nx, m = n)) {

    mat_L <- 1 - tcrossprod(x = data, y = x[ind, , drop = FALSE])
    kde[ind] <- ch * colMeans(L(mat_L / h^2))

  }

  # Result
  return(kde)

}


#' @rdname kde_dir
#' @export
c_h <- function(h, q, L = NULL, log = FALSE) {

  # Von Mises--Fisher kernel?
  if (is.null(L)) {

    # Exponentially scaled von Mises--Fisher normalizing constant:
    # c_h(h, q) = exp(1 / h^2) * c_q(1 / h^2), computed directly, as adding
    # 1 / h^2 to log c_q(1 / h^2) cancels digits for small h.
    log_ch <- log_c_vmf(q = q, kappa = 1 / h^2, scaled = TRUE)

  } else {

    # Numerical integration, with the area of S^{q - 1} coming from the
    # tangent-normal change of variables.
    ch <- sapply(h, function(hh) {

      integrate(f = function(t) L((1 - t) / hh^2) * (1 - t^2)^(q / 2 - 1),
                lower = -1 + 1e-5, upper = 1 - 1e-5, rel.tol = 1e-5,
                subdivisions = 200, stop.on.error = FALSE)$value *
        w_q(q = q - 1)

    })
    log_ch <- -log(ch)

  }

  # Normalizing constant
  return(switch(log + 1, exp(log_ch), log_ch))

}


#' @title Moments of a directional kernel
#'
#' @description Moments \eqn{\lambda_q(L)}, \eqn{b_q(L)}, and \eqn{d_q(L)} of
#' the kernel \code{L} that enter the AMISE bandwidth of \code{bw_dir_ami}.
#'
#' @inheritParams kde_dir
#' @return The moment, a scalar.
#' @noRd
lambda_L <- function(L = NULL, q) {

  # Von Mises--Fisher kernel?
  if (is.null(L)) {

    # Analytical expression
    lambda <- (2 * pi)^(q / 2)

  } else {

    # Numerical integration
    lambda <- 2^(q / 2 - 1) * w_q(q = q - 1) *
      integrate(function(r) L(r) * r^(q / 2 - 1), lower = 0, upper = 100,
                rel.tol = 1e-5, subdivisions = 200, stop.on.error = FALSE)$value

  }

  # Moment
  return(lambda)

}


#' @noRd
b_L <- function(L = NULL, q) {

  # Von Mises--Fisher kernel?
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


#' @noRd
d_L <- function(L = NULL, q) {

  # Von Mises--Fisher kernel?
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


#' @title Leave-one-out kernel density estimates
#'
#' @description Leave-one-out kernel density estimator evaluated at each
#' observation, for a kernel \code{L} and a grid of bandwidths.
#'
#' @inheritParams kde_dir
#' @inheritParams bw_dir_cv
#' @return A matrix of size \code{c(nrow(data), length(h_grid))}.
#' @noRd
loo_kde <- function(data, h_grid, L) {

  # Normalizing constants
  n <- nrow(data)
  ch <- c_h(h = h_grid, q = ncol(data) - 1, L = L)

  # Leave-one-out kernel sums, accumulated over the row blocks with the
  # bandwidths inside, so the memory is O(n * b).
  f_hi <- matrix(nrow = n, ncol = length(h_grid))
  for (ind in row_blocks(n = n)) {

    mat_L <- 1 - tcrossprod(x = data[ind, , drop = FALSE], y = data)
    for (hi in seq_along(h_grid)) {

      L_hi <- (ch[hi] / (n - 1)) * L(mat_L / h_grid[hi]^2)
      L_hi[cbind(seq_along(ind), ind)] <- 0
      f_hi[ind, hi] <- rowSums(L_hi)

    }

  }

  # Result
  return(f_hi)

}



#' @title Integration routines
#'
#' @description Several quadrature rules for integration of functions on
#' \eqn{S^1}, \eqn{S^2}, and \eqn{S^q}, \eqn{q \ge 3}.
#'
#' @param f function to be integrated on \eqn{S^q}. Must be vectorized and
#' accept matrix inputs of size \code{c(nx, q + 1)}.
#' @param N Defaults to \code{5e2}.
#' @param na.rm ignore possible \code{NA}s arising from the evaluation of
#' \code{f}? Defaults to \code{TRUE}.
#' @param f_vect can \code{f} be called in a vectorized form, with matrix
#' input? Defaults to \code{TRUE}.
#' @inheritParams kde_dir
#' @param M number of Monte Carlo replicates. Defaults to \code{1e5}.
#' @param ... further arguments passed to \code{f}.
#' @return A scalar approximating the integral.
#' @details
#' \code{int_cir} is an extension of equation (4.1.11) in Press et al. (1997),
#' a periodic trapezoidal rule. \code{int_sph} employs the
#' \link[=lebedev]{Lebedev quadrature} on \eqn{S^2}. \code{int_hypsph}
#' implements a Monte Carlo integration on \eqn{S^q}.
#' @examples
#' # S^1, trapezoidal rule
#' f <- function(x) rotasym::d_vMF(x = x, mu = c(0, 1), kappa = 2)
#' int_cir(f = f)
#'
#' # S^2, Lebedev rule
#' f <- function(x) rotasym::d_vMF(x = x, mu = c(0, 0, 1), kappa = 2)
#' int_sph(f = f)
#'
#' # S^2, Monte Carlo
#' f <- function(x) rotasym::d_vMF(x = x, mu = c(0, 0, 1), kappa = 2)
#' int_hypsph(f = f, q = 2)
#' @references
#' Lebedev, V. I. and Laikov, D. N. (1999). A quadrature formula for the
#' sphere of the 131st algebraic order of accuracy. \emph{Doklady Mathematics},
#' 59(3):477--481.
#'
#' Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery B. P. (1997).
#' \emph{Numerical Recipes in Fortran 77: The Art of Scientific Computing}.
#' Volume 1. Cambridge University Press, Cambridge. Second edition.
#' @name int


#' @rdname int
#' @export
int_cir <- function(f, N = 5e2, na.rm = TRUE, f_vect = TRUE, ...) {

  # Extension of equation (4.1.11) in Numerical Recipes in F77
  cir_grid <- to_cir(seq(0, 2 * pi, l = N + 1)[-c(N + 1)])
  f_evals <- switch(f_vect + 1, apply(cir_grid, 1, f, ...), f(cir_grid, ...))
  int <- mean(f_evals, na.rm = na.rm) * 2 * pi
  return(int)

}


#' @rdname int
#' @export
int_sph <- function(f, na.rm = TRUE, f_vect = TRUE, ...) {

  # Lebedev quadrature rule
  f_evals <- switch(f_vect + 1, apply(DirStats::lebedev$xyz, 1, f, ...),
                    f(DirStats::lebedev$xyz, ...))
  return(4 * pi * sum(DirStats::lebedev$w * f_evals, na.rm = na.rm))

}


#' @rdname int
#' @export
int_hypsph <- function(f, q, M = 1e5, na.rm = TRUE, f_vect = TRUE, ...) {

  # Monte Carlo approximation
  sims <- rotasym::r_unif_sphere(n = M, p = q + 1)
  f_evals <- switch(f_vect + 1, apply(sims, 1, f, ...), f(sims, ...))
  return(mean(f_evals, na.rm = TRUE) * rotasym::w_p(p = q + 1, log = FALSE))

}


#' @title Convenience functions
#'
#' @description Normalization of data in \eqn{R^{q + 1}} to \eqn{S^q}.
#' Transformations between \eqn{S^1} and \eqn{[0, 2\pi)}, and between
#' \eqn{S^2} and \eqn{[0, 2\pi) \times [0, \pi]}.
#'
#' @param x matrix or vector, in \eqn{S^1} for \code{to_cir}.
#' @param th vector of angles in \eqn{[0, 2\pi)}.
#' @param ph vector of angles in \eqn{[0, \pi]}.
#' @return Euclidean norm (\code{norm}) and normalized data (\code{normalize}).
#' Position in \eqn{S^1} (\code{to_cir}) or in \eqn{[0, 2\pi)} (\code{to_rad}).
#' Position in \eqn{S^2} (\code{to_sph}) or in \eqn{[0, 2\pi) \times [0, \pi]}
#' (\code{to_rad}).
#' @examples
#' # Normalization
#' x <- 1:3
#' norm2(x)
#' normalize(x)
#' x <- rbind(1:3, 3:1)
#' norm2(x)
#' normalize(x)
#'
#' # Circular transformations
#' th <- 1
#' x <- c(0, 1)
#' to_rad(to_cir(th))
#' to_rad(to_cir(c(th, th + 1)))
#' to_cir(to_rad(x))
#' to_cir(to_rad(rbind(x, -x)))
#'
#' # Spherical transformations
#' th <- 2
#' ph <- 1
#' x <- c(0, 1, 0)
#' to_rad(to_sph(th, ph))
#' to_rad(to_sph(c(th, th + 1),
#'               c(ph, ph + 1)))
#' to_sph(to_rad(x)[, 1], to_rad(x)[, 2])
#' to_sph(to_rad(rbind(x, -x))[, 1], to_rad(rbind(x, -x))[, 2])
#' @name conv


#' @rdname conv
#' @export
norm2 <- function(x) {

  if (!is.matrix(x)) x <- t(x)
  return(sqrt(rowSums(x^2)))

}


#' @rdname conv
#' @export
normalize <- function(x) {

  return(x / norm2(x))

}


#' @rdname conv
#' @export
to_cir <- function(th) {

  return(cbind(cos(th), sin(th)))

}


#' @rdname conv
#' @export
to_rad <- function(x) {

  if (!is.matrix(x)) x <- t(x)
  q <- ncol(x) - 1
  if (q > 2) stop("to_rad only works with circular or spherical data")
  rad <- switch(q,
    atan2(x[, 2], x[, 1]) %% (2 * pi),
    cbind(atan2(x[, 2], x[, 1]) %% (2 * pi), acos(x[, 3]) %% pi)
  )
  return(rad)

}


#' @rdname conv
#' @export
to_sph <- function(th, ph) {

  return(cbind(cbind(cos(th), sin(th)) * sin(ph), cos(ph)))

}

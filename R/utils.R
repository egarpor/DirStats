

#' @title Integration routines
#'
#' @description Several quadrature rules for integration of functions on
#' \eqn{\mathbb{S}^1}, \eqn{\mathbb{S}^2}, and \eqn{\mathbb{S}^q},
#' \eqn{q \ge 3}.
#'
#' @param f function to be integrated on \eqn{\mathbb{S}^q}. Must be vectorized
#' and accept matrix inputs of size \code{c(nx, q + 1)}.
#' @param N Defaults to \code{5e2}.
#' @param na.rm ignore possible \code{NA}s arising from the evaluation of
#' \code{f}? Defaults to \code{TRUE}.
#' @param f_vect can \code{f} be called in a vectorized form, with matrix input?
#' Defaults to \code{TRUE}.
#' @inheritParams kde_dir
#' @param M number of Monte Carlo replicates. Defaults to \code{1e5}.
#' @param sim proposal for the Monte Carlo integration in \code{int_hypsph}:
#' either a function of \code{M} returning a matrix of size
#' \code{c(M, q + 1)}, or a matrix of draws on \eqn{\mathbb{S}^q}. Passing a
#' matrix reuses the same sample across several integrations, which is required
#' if the integral is to be a deterministic function of other parameters (as in
#' an optimization). \code{M} is ignored if \code{sim} is a matrix. Defaults to
#' \code{NULL}, which draws uniformly on \eqn{\mathbb{S}^q}.
#' @param f_sim density of \code{sim} with respect to the surface measure, used
#' as the importance weights. Defaults to \code{NULL}, which assumes that
#' \code{sim} is uniform on \eqn{\mathbb{S}^q}.
#' @param ... further arguments passed to \code{f}.
#' @return A scalar approximating the integral.
#' @details
#' \code{int_cir} is an extension of equation (4.1.11) in Press et al. (1997),
#' a periodic trapezoidal rule. \code{int_sph} employs the
#' \link[=lebedev]{Lebedev quadrature} on \eqn{\mathbb{S}^2}. \code{int_hypsph}
#' implements a Monte Carlo integration on \eqn{\mathbb{S}^q}.
#'
#' Uniform draws spread the Monte Carlo effort evenly over
#' \eqn{\mathbb{S}^q}, and therefore miss integrands that concentrate their
#' mass on a small region: the estimate then collapses to (almost) zero. In that
#' situation, supply through \code{sim} a proposal that charges the region where
#' the integrand lives, together with its density \code{f_sim}, and
#' \code{int_hypsph} performs importance sampling instead. Mixing the proposal
#' with a uniform component is advisable, as it bounds the importance weights
#' wherever the proposal is negligible but the integrand is not.
#' @examples
#' # S^1, trapezoidal rule
#' f <- function(x) d_mixvmf(x = x, mu = c(0, 1), kappa = 2, p = 1)
#' int_cir(f = f)
#'
#' # S^2, Lebedev rule
#' f <- function(x) d_mixvmf(x = x, mu = c(0, 0, 1), kappa = 2, p = 1)
#' int_sph(f = f)
#'
#' # S^2, Monte Carlo
#' f <- function(x) d_mixvmf(x = x, mu = c(0, 0, 1), kappa = 2, p = 1)
#' int_hypsph(f = f, q = 2)
#'
#' # S^2, Monte Carlo with a concentrated integrand. Uniform draws underestimate
#' # the integral, while importance sampling from a proposal that mixes the
#' # concentrated density with a uniform component recovers it
#' mu <- rbind(c(0, 0, 1))
#' kappa <- 500
#' f <- function(x) d_mixvmf(x = x, mu = mu, kappa = kappa, p = 1)^2
#' int_hypsph(f = f, q = 2, M = 1e4)
#' M_mix <- 5e3
#' sim <- rbind(movMF::rmovMF(n = M_mix, theta = kappa * mu, alpha = 1),
#'              normalize(matrix(rnorm(n = M_mix * 3), nrow = M_mix, ncol = 3)))
#' f_sim <- function(x) {
#'   0.5 * d_mixvmf(x = x, mu = mu, kappa = kappa, p = 1) + 0.5 / w_q(q = 2)
#' }
#' int_hypsph(f = f, q = 2, sim = sim, f_sim = f_sim)
#' @references
#' Lebedev, V. I. and Laikov, D. N. (1999). A quadrature formula for the
#' sphere of the 131st algebraic order of accuracy. \emph{Doklady Mathematics},
#' 59(3):477--481.
#'
#' Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery B. P. (1997).
#' \emph{Numerical Recipes in Fortran 77: The Art of Scientific Computing}.
#' Volume 1. Cambridge University Press, Cambridge. Second edition.
#' @seealso \code{\link{conv}} for the coordinate transformations used to build
#' the integration grids, and \code{\link{lebedev}} for the quadrature data
#' employed by \code{int_sph}.
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
int_hypsph <- function(f, q, M = 1e5, na.rm = TRUE, f_vect = TRUE, sim = NULL,
                       f_sim = NULL, ...) {

  # Uniform draws, obtained by projecting a N(0, I_{q + 1}) onto S^q, unless a
  # proposal is given
  if (is.null(sim)) {

    sim <- normalize(matrix(rnorm(n = M * (q + 1)), nrow = M, ncol = q + 1))

  } else if (is.function(sim)) {

    sim <- sim(M)

  }

  # Check dimension
  stopifnot(ncol(sim) == q + 1)

  # Monte Carlo approximation
  f_evals <- switch(f_vect + 1, apply(sim, 1, f, ...), f(sim, ...))
  if (is.null(f_sim)) {

    # A uniform proposal has constant density 1 / w_q, which factors out
    return(mean(f_evals, na.rm = na.rm) * w_q(q = q))

  } else {

    # Importance weights
    return(mean(f_evals / f_sim(sim), na.rm = na.rm))

  }

}


#' @title Convenience functions
#'
#' @description Normalization of data in \eqn{\mathbb{R}^{q + 1}} to
#' \eqn{\mathbb{S}^q}. Transformations between \eqn{\mathbb{S}^1} and
#' \eqn{[0, 2\pi)}, and between \eqn{\mathbb{S}^2} and
#' \eqn{[0, 2\pi) \times [0, \pi]}.
#'
#' @param x matrix or vector, in \eqn{\mathbb{S}^1} for \code{to_cir}.
#' @param th vector of angles in \eqn{[0, 2\pi)}.
#' @param ph vector of angles in \eqn{[0, \pi]}.
#' @return
#' \itemize{
#'   \item \code{norm2}: the Euclidean norms of the rows of \code{x}, a vector
#'   of size \code{nrow(x)}.
#'   \item \code{normalize}: the row-normalized \code{x}, a matrix with the
#'   same size as \code{x}.
#'   \item \code{to_cir}: the position in \eqn{\mathbb{S}^1}, a matrix of size
#'   \code{c(length(th), 2)}.
#'   \item \code{to_sph}: the position in \eqn{\mathbb{S}^2}, a matrix of size
#'   \code{c(length(th), 3)}.
#'   \item \code{to_rad}: the angles of \code{x}, a vector of size
#'   \code{nrow(x)} in \eqn{[0, 2\pi)} if \code{x} is in \eqn{\mathbb{S}^1},
#'   or a matrix of size \code{c(nrow(x), 2)} in \eqn{[0, 2\pi) \times
#'   [0, \pi]} if \code{x} is in \eqn{\mathbb{S}^2}.
#' }
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
#' @seealso \code{\link{int}} for integration on \eqn{\mathbb{S}^q}, which
#' relies on these conversions.
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
    # The polar angle acos(x[, 3]) is already in [0, pi]; applying "%% pi"
    # would wrongly wrap the south pole (acos(-1) = pi) back to 0
    cbind(atan2(x[, 2], x[, 1]) %% (2 * pi), acos(x[, 3]))
  )
  return(rad)

}


#' @rdname conv
#' @export
to_sph <- function(th, ph) {

  return(cbind(cbind(cos(th), sin(th)) * sin(ph), cos(ph)))

}


#' @title Surface area of the hypersphere
#'
#' @description Surface area of the hypersphere
#' \eqn{\mathbb{S}^q := \{x \in \mathbb{R}^{q + 1} : ||x|| = 1\}},
#' \eqn{q \ge 0}, given by
#' \deqn{\omega_q := 2\pi^{(q + 1) / 2} / \Gamma((q + 1) / 2).}
#'
#' @param q dimension of \eqn{\mathbb{S}^q}, \eqn{q \ge 0}. A vector.
#' @param log return the logarithm of the area? Defaults to \code{FALSE}.
#' @return The surface area of \eqn{\mathbb{S}^q}, a vector of size
#' \code{length(q)}.
#' @details
#' \eqn{\mathbb{S}^0 = \{-1, 1\}} and its "surface area" is \eqn{2}. The
#' function is vectorized on \code{q}.
#'
#' \code{w_q} is the reciprocal of the density of the uniform distribution on
#' \eqn{\mathbb{S}^q}, and it is the \eqn{q}-parametrization of
#' \code{\link[rotasym]{w_p}}: \code{w_q(q = q)} equals
#' \code{rotasym::w_p(p = q + 1)}.
#' @seealso \code{\link{int_hypsph}} for the Monte Carlo integration that
#' employs the area, and \code{\link{log_c_vmf}} for the von Mises--Fisher
#' normalizing constant that reduces to \eqn{1 / \omega_q} if \eqn{\kappa = 0}.
#' @examples
#' # Areas of S^0, S^1, and S^2
#' w_q(q = 0:2)
#'
#' # Area as a function of q
#' q <- 0:19
#' plot(q, w_q(q = q), type = "o", pch = 16, xlab = "q", ylab = "Area",
#'      main = expression("Surface area of " * S^q))
#'
#' # The area vanishes for large q, hence the log scale
#' w_q(q = 200)
#' w_q(q = 200, log = TRUE)
#' @export
w_q <- function(q, log = FALSE) {

  # Area of S^q: 2 * pi^((q + 1) / 2) / Gamma((q + 1) / 2)
  log_w_q <- log(2) + 0.5 * (q + 1) * log(pi) - lgamma(0.5 * (q + 1))
  return(switch(log + 1, exp(log_w_q), log_w_q))

}


#' @title Defensive mixture proposal for a mixture of von Mises--Fisher
#' densities
#'
#' @description Builds the proposal used to integrate, by importance sampling,
#' integrands whose mass is driven by a mixture of von Mises--Fisher densities.
#'
#' @inheritParams d_mixvmf
#' @inheritParams int
#' @param p_unif proportion of uniform draws in the proposal. Defaults to
#' \code{0.5}.
#' @param kappa_max largest concentration used to build the proposal. Defaults
#' to \code{1e10}.
#' @return A list with the sampled proposal (\code{sim}) and its density
#' (\code{f_sim}), ready to be passed to \code{\link{int_hypsph}}.
#' @noRd
r_defensive_vmf <- function(q, mu, kappa, p, M = 1e4, p_unif = 0.5,
                            kappa_max = 1e10) {

  # Wood's rejection sampler in movMF::rmovMF() does not return for the
  # large-but-finite concentrations that collapsed mixture components produce
  # (the onset is about 1e16 on S^2, growing mildly with the dimension; kappa =
  # Inf errors instead of hanging). Capping is safe, as importance sampling
  # only needs the proposal to cover the target: it widens the proposal and
  # leaves the estimator unbiased, and the integrand keeps the true
  # concentrations. But the cap must clear by a wide margin the concentrations
  # of genuine fits, or the proposal misses the peaks of the target and the
  # integration silently degenerates: a cap of 250, the analogue of the
  # bic_vmf_mix() screening, doubles the selected bandwidth on ordinary data.
  # The default sits six orders of magnitude below the onset of the hang and
  # resolves any bandwidth down to h = 1e-5, as concentrations beyond 1 / h^2
  # are indistinguishable in the integrand
  if (any(kappa > kappa_max)) {

    warning(paste0("Concentrations above kappa_max = ", kappa_max,
                   " were capped to build the proposal; the mixture has a ",
                   "collapsed component."))

  }
  kappa_prop <- pmin(kappa, kappa_max)

  # A p_unif fraction of uniform draws: the mixture part covers the peaks of
  # the integrand, which uniform draws miss entirely when the mixture is
  # concentrated, while the uniform part bounds the importance weights wherever
  # the mixture is negligible but the integrand is not
  M_unif <- round(p_unif * M)
  M_mix <- M - M_unif
  sim <- rbind(movMF::rmovMF(n = M_mix, theta = kappa_prop * mu, alpha = p),
               normalize(matrix(rnorm(n = M_unif * (q + 1)), nrow = M_unif,
                                ncol = q + 1)))

  # Density of the proposal, with the realized proportions of the two parts and
  # the same capped concentrations as the sampler: otherwise the importance
  # weights do not correspond to the draws
  w <- w_q(q = q)
  f_sim <- function(x) {

    (M_mix / M) * d_mixvmf(x = x, mu = mu, kappa = kappa_prop, p = p) +
      (M_unif / M) / w

  }

  # Result
  return(list(sim = sim, f_sim = f_sim))

}



#' @title Cross-validation bandwidth selectors for directional data
#'
#' @description Likelihood and least squares cross-validation bandwidth
#' selectors for kernel density estimation with directional data.
#'
#' @param h_grid vector of bandwidths for performing a grid search. Defaults
#' to\cr \code{exp(seq(log(0.05), log(1.5), l = 100))}.
#' @inheritParams kde_dir
#' @param plot_it display an informative plot on the optimization's grid search?
#' Defaults to \code{FALSE}.
#' @param optim run an optimization? Defaults to \code{TRUE}. Otherwise, a grid
#' search on \code{h} is done. For \code{bw_dir_lcv} and \code{bw_dir_lscv},
#' only effective if \code{L = NULL}.
#' @param optim_par,optim_lower,optim_upper initial value, lower bound, and
#' upper bound of the bandwidth for \code{\link[stats]{optim}} when using the
#' \code{"L-BFGS-B"} method. The default \code{optim_lower} of
#' \code{bw_dir_lscv} and \code{bw_dir_emi} prevents the overflow of their
#' losses for large \eqn{q}.
#' @return A list with entries:
#' \itemize{
#'   \item \code{h_opt}: selected bandwidth.
#'   \item \code{h_grid}: \code{h_grid}, if used (otherwise \code{NULL}).
#'   \item \code{CV_opt}: minimum of the CV loss.
#'   \item \code{CV_grid}: value of the CV function at \code{h_grid}, if used
#'   (otherwise \code{NULL}).
#' }
#' @details
#' \code{data} is not checked to have unit norm, so the user must be careful.
#'
#' When \code{L = NULL}, faster vectorized expressions specific to the von
#' Mises--Fisher kernel are employed.
#'
#' The pairwise sums are accumulated in row blocks of at most 2.5e5 inner
#' products, so the memory footprint stays flat in \eqn{n} for \eqn{n > 5000}.
#'
#' \code{bw_dir_lscv} employs Monte Carlo integration for \eqn{q > 2}, which
#' results in a random output. Use \code{set.seed} before to avoid it.
#'
#' A warning is issued if \code{optim} does not converge, or if the selected
#' bandwidth is at an extreme of \code{h_grid} or of
#' \code{c(optim_lower, optim_upper)}.
#' @section Warning:
#' \code{bw_dir_lscv} is not suited for large dimensions, as its loss overflows
#' for small bandwidths from \eqn{q \approx 45}. Use instead
#' \href{https://CRAN.R-project.org/package=polykde}{\code{polykde::bw_cv_polysph}}
#' with \code{type = "LSCV"} and \code{arcsinh = TRUE}.
#' @references
#' Hall, P., Watson, G. S., and Cabrera, J. (1987). Kernel density estimation
#' with spherical data. \emph{Biometrika}, 74(4):751--762.
#' \doi{10.1093/biomet/74.4.751}
#' @seealso \code{\link{kde_dir}} for the kernel density estimator these
#' bandwidths are intended for, and \code{\link{bw_dir_pi}} for plug-in
#' alternatives.
#' @examples
#' # Sample
#' n <- 25
#' q <- 2
#' if (requireNamespace("rotasym", quietly = TRUE)) {
#'
#'   set.seed(42)
#'   samp <- rotasym::r_vMF(n = n, mu = c(1, rep(0, q)), kappa = 2)
#'
#'   # bw_dir_lcv
#'   print(bw_dir_lcv(data = samp, optim = TRUE)$h_opt)
#'   print(bw_dir_lcv(data = samp, optim = FALSE, plot_it = TRUE)$h_opt)
#'   print(bw_dir_lcv(data = samp, L = function(x) exp(-x))$h_opt)
#'
#'   # bw_dir_lscv
#'   set.seed(42)
#'   print(bw_dir_lscv(data = samp, optim = TRUE)$h_opt)
#'   print(bw_dir_lscv(data = samp, optim = FALSE, plot_it = TRUE)$h_opt)
#'   print(bw_dir_lscv(data = samp, L = function(x) exp(-x))$h_opt)
#'
#' }
#' @name bw_dir_cv


#' @rdname bw_dir_cv
#' @export
bw_dir_lcv <- function(data, h_grid = exp(seq(log(0.05), log(1.5), l = 100)),
                       L = NULL, plot_it = FALSE, optim = TRUE,
                       optim_par = 0.25, optim_lower = 1e-4, optim_upper = 10) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Dimensions of data
  n <- nrow(data)
  q <- ncol(data) - 1

  # Row blocks holding up to 2.5e5 inner products each, so the pairwise sums
  # take O(n * b) memory instead of the O(n^2) of all the pairs at once.
  blocks <- row_blocks(n = n)

  # Type of approach
  if (is.null(L)) {

    # Row-wise log-sum-exp shifts m[i] = max_{j != i} (X_i'X_j - 1), computed
    # once outside the bandwidth loop, as those of (X_i'X_j - 1) / h^2 are just
    # m / h^2.
    m <- numeric(n)
    for (ind in blocks) {

      M_b <- tcrossprod(x = data[ind, , drop = FALSE], y = data) - 1
      M_b[cbind(seq_along(ind), ind)] <- -Inf
      m[ind] <- M_b[cbind(seq_along(ind), max.col(M_b, ties.method = "first"))]

    }

    # Shifted leave-one-out log-kernel matrix of a row block:
    # log exp((X_i'X_j - 1) / h^2) = (M_b[i, j] + m[i]) / h^2, with -Inf on the
    # diagonal. Recomputing it per evaluation costs 1.45 times the exp() it
    # feeds, so it is cached whole for n <= 5000 (200 Mb) and recomputed in
    # blocks above, keeping the memory flat in n.
    M_block <- function(ind) {

      M_b <- tcrossprod(x = data[ind, , drop = FALSE], y = data) - (1 + m[ind])
      M_b[cbind(seq_along(ind), ind)] <- -Inf
      return(M_b)

    }
    if (n <= 5000) {

      blocks <- list(seq_len(n))
      M <- M_block(ind = blocks[[1]])
      M_block <- function(ind) M

    }

    # CV function, reusing the blocks and M
    CV_vmf <- function(h) {

      ch <- c_h(h = h, q = q, log = TRUE) - log(n - 1)
      cv <- numeric(length(h))
      for (ind in blocks) {

        M_b <- M_block(ind = ind)
        m_b <- m[ind]
        cv <- cv + sapply(seq_along(h), function(hi) {

          h2 <- h[hi]^2
          sum(ch[hi] + m_b / h2 + log(rowSums(exp(M_b / h2))))

        })

      }
      -cv / n

    }

    loss <- CV_vmf

  } else {

    # Loss on h_grid, which is always searched for a custom kernel
    CV <- -colMeans(log(loo_kde(data = data, h_grid = h_grid, L = L)))
    loss <- function(h) CV
    optim <- FALSE

  }

  # Result
  return(bw_select(loss = loss, optim = optim, h_grid = h_grid,
                   optim_par = optim_par, optim_lower = optim_lower,
                   optim_upper = optim_upper, plot_it = plot_it,
                   ylab = "LCV", name = "CV"))

}


#' @rdname bw_dir_cv
#' @export
bw_dir_lscv <- function(data, h_grid = exp(seq(log(0.05), log(1.5), l = 100)),
                        L = NULL, plot_it = FALSE, optim = TRUE,
                        optim_par = 0.25,
                        optim_lower = max(1e-4, exp(-700 / (ncol(data) - 1)) /
                                            sqrt(2 * pi)),
                        optim_upper = 10) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Dimensions of data
  n <- nrow(data)
  q <- ncol(data) - 1

  # Row blocks holding up to 2.5e5 inner products each, so the pairwise sums
  # take O(n * b) memory instead of the O(n^2) of all the pairs at once.
  blocks <- row_blocks(n = n)

  # Type of approach
  if (is.null(L)) {

    # Pairs (i, j) with j < i of the b-th row block, the strict lower triangle
    # held by the columns up to max(ind): X_i'X_j - 1 and ||X_i + X_j||. They
    # are cached for n <= 5000, keeping the memory flat in n.
    pairs_block <- function(b) {

      ind <- blocks[[b]]
      j_max <- ind[length(ind)]
      d_ij <- tcrossprod(x = data[ind, , drop = FALSE],
                         y = data[seq_len(j_max), , drop = FALSE])
      d_ij <- d_ij[outer(X = ind, Y = seq_len(j_max), FUN = ">")]
      return(list(d_ij = d_ij - 1, norm_ij = sqrt(pmax(2 + 2 * d_ij, 0))))

    }
    if (n <= 5000) {

      pairs_cache <- lapply(seq_along(blocks), pairs_block)
      pairs_block <- function(b) pairs_cache[[b]]

    }

    # CV function, reusing the blocks and their pairs
    CV_vmf <- function(h) {

      # c_h() is by construction the exponentially scaled vMF constant,
      # c_h(h, q) = exp(1 / h^2) * c_q(1 / h^2) =: \tilde c_q(1 / h^2), so
      # log_Cq = log \tilde c_q(1 / h^2) and log_Cq2 = log \tilde c_q(2 / h^2).
      # Their exp() factors cancel the exponentially small terms in the sums
      # below. Everything stays in log scale until each complete term is
      # formed, as Cq^2 overflows for large q and small h.
      log_Cq <- c_h(h = h, q = q, log = TRUE)
      log_Cq2 <- c_h(h = h / sqrt(2), q = q, log = TRUE)

      # Accumulate the logs of the pairwise kernel sums block by block, with
      # exponents that never overflow:
      # (i) exp((X_i'X_j - 1) / h^2) <= 1, as X_i'X_j <= 1;
      # (ii) with kappa_ij = ||X_i + X_j|| / h^2 and
      # c_q(kappa) = exp(-kappa) * \tilde c_q(kappa),
      #  exp(-2 / h^2) / c_q(kappa_ij)
      #   = exp((||X_i + X_j|| - 2) / h^2) / \tilde c_q(kappa_ij),
      # whose exponent is <= 0 as ||X_i + X_j|| <= 2. log_c_vmf(scaled = TRUE)
      # gives log \tilde c_q(kappa_ij) directly, so kappa_ij is never added.
      log_S1 <- log_S2 <- matrix(nrow = length(blocks), ncol = length(h))
      for (b in seq_along(blocks)) {

        pairs_b <- pairs_block(b = b)
        log_S1[b, ] <- sapply(h, function(hh) {

          log_sum_exp(pairs_b$d_ij / hh^2)

        })
        log_S2[b, ] <- sapply(h, function(hh) {

          log_sum_exp((pairs_b$norm_ij - 2) / hh^2 -
                        log_c_vmf(q = q, kappa = pairs_b$norm_ij / hh^2,
                                  scaled = TRUE))

        })

      }
      log_S1 <- apply(log_S1, 2, log_sum_exp)
      log_S2 <- apply(log_S2, 2, log_sum_exp)

      # Computation of the CV, exponentiating each complete term
      exp(2 * log_Cq - log_Cq2 - log(n)) +
        exp(log(2) + 2 * log_Cq - 2 * log(n) + log_S2) -
        exp(log(4) + log_Cq - log(n) - log(n - 1) + log_S1)

    }

    loss <- CV_vmf

  } else {

    # \hat f_(h, -i)(X_i) terms
    f_hi <- loo_kde(data = data, h_grid = h_grid, L = L)

    # Integral term. For q > 2, importance sampling from the estimator itself,
    # a mixture of von Mises--Fisher densities centered at the data, as uniform
    # draws miss its mass for small bandwidths. The proposal is a defensive
    # mixture of the estimator with a uniform.
    int <- sapply(seq_along(h_grid), function(hi) {

      prop <- switch((q > 2) + 1, NULL,
                     safe_mixvmf(q = q, mu = data,
                                 kappa = rep(1 / h_grid[hi]^2, n),
                                 p = rep(1 / n, n), M = 1e4))
      int_dir(f = function(x) kde_dir(x = x, data = data, h = h_grid[hi],
                                      L = L)^2,
              q = q, sim = prop$sim, f_sim = prop$f_sim)

    })

    # Loss on h_grid, which is always searched for a custom kernel
    CV <- int - 2 / n * colSums(f_hi)
    loss <- function(h) CV
    optim <- FALSE

  }

  # Result
  return(bw_select(loss = loss, optim = optim, h_grid = h_grid,
                   optim_par = optim_par, optim_lower = optim_lower,
                   optim_upper = optim_upper, plot_it = plot_it,
                   ylab = "LSCV", name = "CV"))

}



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
#' search on \code{h} is done. Only effective if \code{L = NULL}.
#' @param optim_par,optim_lower,optim_upper parameters passed to \code{par},
#' \code{lower}, and \code{upper} in \code{\link[stats]{optim}} when using the
#' \code{"L-BFGS-B"} method. Default to \code{0.25},
#' \code{sqrt(.Machine$double.eps)}, and \code{10}.
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
#' When \code{L = NULL}, faster vectorized expressions specific to the von
#' Mises--Fisher kernel are employed. The pairwise sums are accumulated in row
#' blocks of at most 2.5e5 inner products, so the memory footprint stays flat
#' in \eqn{n}, except that \code{bw_dir_lcv} caches the whole \eqn{n \times n}
#' matrix for \eqn{n \le 5000} for speed.
#'
#' \code{bw_dir_lscv} employs Monte Carlo integration for \eqn{q > 2}, which
#' results in a random output. Use \code{set.seed} before to avoid it.
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
                       optim_par = 0.25, optim_lower =
                         sqrt(.Machine$double.eps), optim_upper = 10) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Dimensions of data
  n <- nrow(data)
  q <- ncol(data) - 1

  # Row blocks holding up to 2.5e5 inner products each, so the pairwise sums
  # take O(n * b) memory instead of the O(n^2) of all the pairs at once
  b <- max(1, min(n, floor(2.5e5 / n)))
  blocks <- split(x = seq_len(n), f = ceiling(seq_len(n) / b))

  # Type of approach
  if (is.null(L)) {

    # Row-wise log-sum-exp shifts m[i] = max_{j != i} (X_i'X_j - 1), hoisted
    # out of the bandwidth loop as those of (X_i'X_j - 1) / h^2 are just
    # m / h^2. Without them the leave-one-out sums underflow to zero for small h
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
    # blocks above, keeping the memory flat in n
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

      ch <- log(c_h(h = h, q = q) / (n - 1))
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

    # Optimize versus grid search
    if (optim) {

      # Optimization of CV
      opt <- optim(fn = CV_vmf, par = optim_par, lower = optim_lower,
                   upper = optim_upper, method = "L-BFGS-B")

      # Result
      return(list(h_opt = opt$par, h_grid = NULL, CV_opt = opt$value,
                  CV_grid = NULL))

    } else {

      CV <- CV_vmf(h = h_grid)

    }

  } else {

    # Normalizing constant
    ch <- c_h(h = h_grid, q = q, L = L)

    # Leave-one-out kernel sums, accumulated over the row blocks with the
    # bandwidths inside, so the memory is O(n * b) instead of the
    # O(n^2 * length(h_grid)) of the full array
    CV <- matrix(nrow = n, ncol = length(h_grid))
    for (ind in blocks) {

      mat_L <- 1 - tcrossprod(x = data[ind, , drop = FALSE], y = data)
      for (hi in seq_along(h_grid)) {

        L_hi <- (ch[hi] / (n - 1)) * L(mat_L / h_grid[hi]^2)
        L_hi[cbind(seq_along(ind), ind)] <- 0
        CV[ind, hi] <- rowSums(L_hi)

      }

    }

    # Computation of the CV
    CV <- -colMeans(log(CV))

  }

  # Index of the minimun CV
  ind_min <- which.min(CV)

  # Plot?
  if (plot_it) {

    plot(h_grid, CV, xlab = "h", ylab = "LCV", type = "o")
    rug(h_grid)
    abline(v = h_grid[ind_min], col = 2)

  }

  # Warning
  if (h_grid[ind_min] == max(h_grid) || h_grid[ind_min] == min(h_grid)) {

    warning("h_opt is at the extreme of h_grid")

  }

  # Result
  return(list(h_opt = h_grid[ind_min], h_grid = h_grid,
              CV_opt = CV[ind_min], CV_grid = CV))

}


#' @rdname bw_dir_cv
#' @export
bw_dir_lscv <- function(data, h_grid = exp(seq(log(0.05), log(1.5), l = 100)),
                        L = NULL, plot_it = FALSE, optim = TRUE,
                        optim_par = 0.25, optim_lower =
                          sqrt(.Machine$double.eps), optim_upper = 10) {

  # Remove NA's
  stopifnot(is.matrix(data))
  data <- na.omit(data)

  # Dimensions of data
  n <- nrow(data)
  q <- ncol(data) - 1

  # Row blocks holding up to 2.5e5 inner products each, so the pairwise sums
  # take O(n * b) memory instead of the O(n^2) of all the pairs at once. The
  # evaluation time is flat in the block size above 1e5, while the peak memory
  # grows linearly with it, as log_c_vmf() allocates about ten temporaries of
  # the block length
  b <- max(1, min(n, floor(2.5e5 / n)))
  blocks <- split(x = seq_len(n), f = ceiling(seq_len(n) / b))

  # Type of approach
  if (is.null(L)) {

    # CV function, reusing the block partition above
    CV_vmf <- function(h) {

      # c_h() is by construction the exponentially scaled vMF constant,
      # c_h(h, q) = exp(1 / h^2) * c_q(1 / h^2) =: \tilde c_q(1 / h^2), so
      # Cq = \tilde c_q(1 / h^2) and Cq2 = \tilde c_q(2 / h^2). Their exp()
      # factors cancel the exponentially small terms in the sums below
      Cq <- c_h(h = h, q = q)
      Cq2 <- c_h(h = h / sqrt(2), q = q)

      # Accumulate the pairwise kernel sums block by block, with exponents
      # that never overflow:
      # (i) exp((X_i'X_j - 1) / h^2) <= 1, as X_i'X_j <= 1;
      # (ii) with kappa_ij = ||X_i + X_j|| / h^2 and
      # c_q(kappa) = exp(-kappa) * \tilde c_q(kappa),
      #  exp(-2 / h^2) / c_q(kappa_ij)
      #   = exp((||X_i + X_j|| - 2) / h^2) / \tilde c_q(kappa_ij),
      # whose exponent is <= 0 as ||X_i + X_j|| <= 2. The scaling is never
      # formed: log_c_vmf() gives log c_q(kappa_ij) and the exp(kappa_ij)
      # stays in the exponent
      sum_CV1 <- sum_CV2 <- numeric(length(h))
      for (ind in blocks) {

        # Strict lower triangle of the row block: the pairs (i, j) with j < i,
        # held by the columns up to max(ind)
        j_max <- ind[length(ind)]
        d_ij <- tcrossprod(x = data[ind, , drop = FALSE],
                           y = data[seq_len(j_max), , drop = FALSE])
        d_ij <- d_ij[outer(X = ind, Y = seq_len(j_max), FUN = ">")]

        # ||X_i + X_j|| = sqrt(2 * (1 + X_i'X_j)), truncated as X_i'X_j may
        # fall below -1 for numerically antipodal observations
        norm_ij <- sqrt(pmax(2 + 2 * d_ij, 0))
        d_ij <- d_ij - 1
        sum_CV1 <- sum_CV1 + sapply(h, function(hh) sum(exp(d_ij / hh^2)))
        sum_CV2 <- sum_CV2 + sapply(h, function(hh) {
          sum(exp(-2 / hh^2 - log_c_vmf(q = q, kappa = norm_ij / hh^2)))
        })

      }

      # Computation of the CV
      -((4 * Cq) / (n * (n - 1)) * sum_CV1 - Cq^2 / (n * Cq2) -
          (2 * Cq^2) / n^2 * sum_CV2)

    }

    # Optimize versus grid search
    if (optim) {

      # Optimization
      opt <- optim(fn = CV_vmf, par = optim_par, lower = optim_lower,
                   upper = optim_upper, method = "L-BFGS-B")

      # Result
      return(list(h_opt = opt$par, h_grid = NULL, CV_opt = opt$value,
                  CV_grid = NULL))

    } else {

      CV <- CV_vmf(h = h_grid)

    }

  } else {

    # Normalizing constant
    ch <- c_h(h = h_grid, q = q, L = L)

    # \hat f_(h, -i)(X_i) terms, accumulated over the row blocks with the
    # bandwidths inside, so the memory is O(n * b) instead of the
    # O(n^2 * length(h_grid)) of the full array
    f_hi <- matrix(nrow = n, ncol = length(h_grid))
    for (ind in blocks) {

      mat_L <- 1 - tcrossprod(x = data[ind, , drop = FALSE], y = data)
      for (hi in seq_along(h_grid)) {

        L_hi <- (ch[hi] / (n - 1)) * L(mat_L / h_grid[hi]^2)
        L_hi[cbind(seq_along(ind), ind)] <- 0
        f_hi[ind, hi] <- rowSums(L_hi)

      }

    }

    # Integral term
    if (q == 1) {

      int <- sapply(seq_along(h_grid), function(hi) {
        int_cir(f = function(x) kde_dir(x = x, data = data, h = h_grid[hi],
                                        L = L)^2,
                N = 500)
      })

    } else if (q == 2) {

      int <- sapply(seq_along(h_grid), function(hi) {
        int_sph(f = function(x) kde_dir(x = x, data = data, h = h_grid[hi],
                                        L = L)^2)
      })

    } else {

      # The estimator is itself a mixture of von Mises--Fisher densities
      # centered at the data, so it gives the importance sampling proposal.
      # Uniform draws miss its mass for small bandwidths, precisely where the
      # criterion must be accurate.
      # ponytail: the proposal is rebuilt for each h, as its concentration is
      # 1 / h^2, so the criterion is slightly jagged along h_grid. Fine for a
      # grid search; use common random numbers if it is ever optimized.
      int <- sapply(seq_along(h_grid), function(hi) {
        prop <- r_defensive_vmf(q = q, mu = data,
                                kappa = rep(1 / h_grid[hi]^2, n),
                                p = rep(1 / n, n), M = 1e4)
        int_hypsph(f = function(x) kde_dir(x = x, data = data, h = h_grid[hi],
                                           L = L)^2,
                   q = q, sim = prop$sim, f_sim = prop$f_sim)
      })

    }

    # CV
    CV <- 2 / n * colSums(f_hi) - int
    CV <- -CV

  }

  # Index of the minimun CV
  ind_min <- which.min(CV)

  # Plot?
  if (plot_it) {

    plot(h_grid, CV, xlab = "h", ylab = "LSCV", type = "o")
    rug(h_grid)
    abline(v = h_grid[ind_min], col = 2)

  }

  # Warning
  if (h_grid[ind_min] == max(h_grid) || h_grid[ind_min] == min(h_grid)) {

    warning("h_opt is at the extreme of h_grid")

  }

  # Result
  return(list(h_opt = h_grid[ind_min], h_grid = h_grid,
              CV_opt = CV[ind_min], CV_grid = CV))

}

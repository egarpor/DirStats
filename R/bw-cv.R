

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
#' Mises--Fisher kernel are employed.
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
#' set.seed(42)
#' samp <- rotasym::r_vMF(n = n, mu = c(1, rep(0, q)), kappa = 2)
#'
#' # bw_dir_lcv
#' bw_dir_lcv(data = samp, optim = TRUE)$h_opt
#' bw_dir_lcv(data = samp, optim = FALSE, plot_it = TRUE)$h_opt
#' bw_dir_lcv(data = samp, L = function(x) exp(-x))$h_opt
#'
#' # bw_dir_lscv
#' set.seed(42)
#' bw_dir_lscv(data = samp, optim = TRUE)$h_opt
#' bw_dir_lscv(data = samp, optim = FALSE, plot_it = TRUE)$h_opt
#' bw_dir_lscv(data = samp, L = function(x) exp(-x))$h_opt
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

  # Type of approach
  if (is.null(L)) {

    # Ensure data types for lcv_dir_vmf
    data_dir <- matrix(as.double(data), nrow = n, ncol = q + 1)
    n <- as.integer(n)
    q <- as.integer(q)

    # Optimize versus grid search
    if (optim) {

      # CV function calling lcv_dir_vmf
      CV <- function(h) {
        .Fortran("lcv_dir_vmf", data_dir = data_dir, h = as.double(h),
                 ch = as.double(c_h(h = h, q = q, L = L)), n = n,
                 nq = q, lenh = as.integer(length(h)),
                 CV = double(length(h)))$CV
      }

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

    # Leave-one-out kernel sums, one bandwidth at a time. This avoids
    # allocating the full n x n x length(h_grid) array (and the subsequent
    # apply()), keeping the memory footprint at O(n^2) instead of
    # O(n^2 * length(h_grid))
    mat_L <- 1 - data %*% t(data)
    CV <- matrix(nrow = n, ncol = length(h_grid))
    for (hi in seq_along(h_grid)) {

      L_hi <- (ch[hi] / (n - 1)) * L(mat_L / h_grid[hi]^2)
      diag(L_hi) <- 0
      CV[, hi] <- rowSums(L_hi)

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

  # Type of approach
  if (is.null(L)) {

    if (!R_code) {

      # Ensure data types for lscv_dir_vmf
      data_dir <- matrix(as.double(data), nrow = n, ncol = p)
      n <- as.integer(n)
      q <- as.integer(q)
      km1 <- as.integer(floor((q - 1) / 2) + 1)
      B <- as.double(rep(0, km1))

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

    # \hat f_(h, -i)(X_i) terms
    f_hi <- matrix(nrow = n, ncol = length(h_grid))
    for (i in 1:n) {
      for (hi in seq_along(h_grid)) {

        f_hi[i, hi] <- kde_dir(x = data[i, ], data = data[-i, ], h = h_grid[hi])

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

      int <- sapply(seq_along(h_grid), function(hi) {
        int_hypsph(f = function(x) kde_dir(x = x, data = data, h = h_grid[hi],
                                           L = L)^2,
                   q = q, M = 1e4)
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

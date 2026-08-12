
# Log-Bessel spline tables consumed by log_besselI_scaled(..., spline = TRUE).
#
# The (scaled, log) Bessel function log(exp(-x) * I_nu(x)) is smooth in
# u = log(x), so it is interpolated on a uniform grid in u. Below x_lo an
# analytic small-argument series is used and above x_hi a large-argument
# asymptotic (Bessel::besselIasym) is used, so the stored spline only needs to
# cover [x_lo, x_hi]. For nu above max(bessel_nus) a large-order asymptotic
# (Bessel::besselI.nuAsym) takes over, so the table stops at that order.
#
# Run with the working directory set to data-raw/ to (re)generate
# ../R/sysdata.rda.

# Grid parameters. The order grid stops at nu = 50: beyond that,
# Bessel::besselI.nuAsym(k.max = 5) is accurate to < 1e-10 and takes over,
# which halves the table relative to going up to nu = 100.
bessel_nus <- seq(0, 50, by = 0.5)
bessel_x_lo <- 1e-5
bessel_x_hi <- 1e4
bessel_N <- 1400

# Uniform knots in u = log(x)
u_bessel <- seq(log(bessel_x_lo), log(bessel_x_hi), length.out = bessel_N)
x_bessel <- exp(u_bessel)

# Analytic leading term of log(exp(-x) * I_nu(x)) as x -> 0
log_besselI_series <- function(x, nu) {
  if (nu == 0) -x else nu * log(x / 2) - lgamma(nu + 1) - x
}

# Evaluate the (scaled, log) Bessel on the grid, one column per order
logI_bessel <- matrix(NA_real_, nrow = bessel_N, ncol = length(bessel_nus))
for (j in seq_along(bessel_nus)) {

  nu <- bessel_nus[j]
  cat("nu =", nu, "\n")
  y <- log(besselI(x = x_bessel, nu = nu, expon.scaled = TRUE))

  # Fill any underflow (-Inf) with the analytic series for small x. This does
  # trigger at the largest orders (e.g. nu = 50 near x = x_lo), where the scaled
  # Bessel underflows to 0; the leading small-argument series is exact there.
  bad <- !is.finite(y)
  if (any(bad)) {

    y[bad] <- log_besselI_series(x_bessel[bad], nu)

  }
  stopifnot(all(is.finite(y)))
  logI_bessel[, j] <- y

}

# Save tables
save(list = c("u_bessel", "bessel_nus", "bessel_x_lo", "bessel_x_hi",
              "logI_bessel"),
     file = "../R/sysdata.rda", compress = "xz")


.onLoad <- function(libname = find.package("DirStats"), pkgname = "DirStats") {

  # Build the (scaled, log) Bessel splines in u = log(x) from the tables
  # (u_bessel, bessel_nus, logI_bessel) stored in sysdata.rda, one splinefun
  # per order, listed so that the order bessel_nus[j] is at position j. These
  # are consumed by log_besselI_scaled(..., spline = TRUE).
  spl <- lapply(seq_along(bessel_nus), function(j) {

    splinefun(x = u_bessel, y = logI_bessel[, j], method = "fmm")

  })
  assign(x = "bessel_splines", value = spl, pos = asNamespace(pkgname))

  # CRAN NOTE avoidance for the sysdata symbols used above
  if (getRversion() >= "2.15.1") {

    utils::globalVariables(c("u_bessel", "bessel_nus", "bessel_x_lo",
                             "bessel_x_hi", "logI_bessel", "bessel_splines"))

  }
  invisible()

}

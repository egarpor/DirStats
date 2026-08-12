
.onLoad <- function(libname = find.package("DirStats"), pkgname = "DirStats") {

  # Build the (scaled, log) Bessel splines in u = log(x) from the tables
  # (u_bessel, bessel_nus, logI_bessel) stored in sysdata.rda, one splinefun
  # per order. These are consumed by log_besselI_scaled(..., spline = TRUE).
  env <- asNamespace(pkgname)
  for (j in seq_along(bessel_nus)) {

    key <- sprintf("%03d", round(10 * bessel_nus[j]))
    assign(x = paste0("log_besselI_scaled_spline_", key),
           value = splinefun(x = u_bessel, y = logI_bessel[, j],
                             method = "fmm"),
           pos = env)

  }

  # CRAN NOTE avoidance for the sysdata symbols used above
  if (getRversion() >= "2.15.1") {

    utils::globalVariables(c("u_bessel", "bessel_nus", "bessel_x_lo",
                             "bessel_x_hi", "logI_bessel"))

  }
  invisible()

}

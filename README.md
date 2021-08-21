# DirStats

[![License:
GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![](https://app.travis-ci.com/egarpor/DirStats.svg?branch=master)](https://app.travis-ci.com/egarpor/DirStats)
[![](https://www.r-pkg.org/badges/version/DirStats?color=green)](https://cran.r-project.org/package=DirStats)
[![](http://cranlogs.r-pkg.org/badges/grand-total/DirStats?color=green)](https://cran.r-project.org/package=DirStats)
[![](http://cranlogs.r-pkg.org/badges/last-month/DirStats?color=green)](https://cran.r-project.org/package=DirStats)

<!-- <img src="" alt="DirStats  hexlogo" align="right" width="200" style="padding: 0 15px; float: right;"/> -->

## Overview

Currently implementing nonparametric kernel density estimation,
bandwidth selection, and other utilities for analyzing directional data.
Further nonparametric tools expected to be included in subsequent
releases.

## Installation

Get the released version from CRAN:

``` r
# Install the package
install.packages("DirStats")

# Load package
library(DirStats)
```

Alternatively, get the latest version from GitHub:

``` r
# Install the package
library(devtools)
install_github("egarpor/DirStats")

# Load package
library(DirStats)
```

## Usage

The following are examples of the usage of the Bai et al. (1988)’s
kernel density estimator, the cross-validatory bandwidth selectors in
Hall et al. (1987), and the plug-in bandwidth selectors in
García-Portugués (2013).

### Compute bandwidths

``` r
# Sample from a von Mises--Fisher on S^2
q <- 2
n <- 300
set.seed(42)
samp <- rbind(rotasym::r_vMF(n = n / 3, mu = c(rep(0, q), 1), kappa = 5),
              rotasym::r_vMF(n = n / 3, mu = c(rep(0, q), -1), kappa = 5),
              rotasym::r_vMF(n = n / 3, mu = c(1, rep(0, q)), kappa = 5))

# LCV bandwidth
(h_LCV <- bw_dir_lcv(data = samp)$h_opt)
#> [1] 0.2567438

# LSCV bandwidth
(h_LSCV <- bw_dir_lscv(data = samp)$h_opt)
#> [1] 0.232965

# ROT bandwidth
(h_ROT <- bw_dir_rot(data = samp))
#> [1] 0.4053648

# Mixture fit, for AMI and EMI bandwidth selectors
fit_mix <- bic_vmf_mix(data = samp)

# AMI bandwidth
(h_AMI <- bw_dir_ami(data = samp, fit_mix = fit_mix))
#> [1] 0.2054242

# EMI bandwidth
(h_EMI <- bw_dir_emi(data = samp, fit_mix = fit_mix)$h_opt)
#> [1] 0.22527
```

### Compute kernel density estimator

``` r
# Compute kde
l <- 200
th <- seq(0, 2 * pi, l = l)
phi <- seq(0, pi, l = l / 2)
ang <- expand.grid(th = th, phi = phi)
x <- to_sph(th = ang$th, ph = ang$phi)
kde <- kde_dir(x = x, data = samp, h = h_EMI)

# Visualization
contour(x = th, y = phi, z = matrix(kde, nrow = l, ncol = l / 2),
        col = viridisLite::viridis(15), 
        levels = round(seq(min(kde), max(kde), length.out = 15), 4), 
        lwd = 2, xlab = expression(theta), ylab = expression(phi))
points(to_rad(samp), pch = 16)
```

<img src="README/README-kde-1.png" style="display: block; margin: auto;" />

## References

Bai, Z. D., Rao, C. R., and Zhao, L. C. (1988). Kernel estimators of
density function of directional data. *Journal of Multivariate
Analysis*, 27(1):24–39.
[doi:10.1016/0047-259X(88)90113-3](https://doi.org/10.1016/0047-259X(88)90113-3).

García-Portugués, E. (2013). Exact risk improvement of bandwidth
selectors for kernel density estimation with directional data.
*Electronic Journal of Statistics*, 7:1655–1685.
[doi:10.1214/13-ejs821](https://doi.org/10.1214/13-ejs821).

Hall, P., Watson, G. S., and Cabrera, J. (1987). Kernel density
estimation with spherical data. *Biometrika*, 74(4):751–762.
[doi:10.1093/biomet/74.4.751](https://doi.org/10.1093/biomet/74.4.751).

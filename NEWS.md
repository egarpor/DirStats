# DirStats 0.1.5

* Initial version.

# DirStats 0.1.6

* Fixed installation problem in r-patched-solaris-x86.

# DirStats 0.1.7

* Fixed bug when calling `int_cir()` in `bw_dir_emi()` (reported by Yong Wang).

# DirStats 0.1.8

* Fixed bug in `bw_rot() `for q = 2 and improve stability in `bw_rot()`.

# DirStats 0.1.9

* Expose `min_kappa` and `max_kappa` in `kappa_ml()`.

# DirStats 0.1.10

* Fix 'Fortran 2018 deleted feature' warning.

# DirStats 1.0.0

* Add `log_besselI_scaled()` and `log_c_vmf()` for the stable evaluation of the scaled log-Bessel function and the log von Mises--Fisher normalizing constant.
* Add `Bessel` to `Imports` and move `rotasym` to `Suggests`.
* Remove all FORTRAN code in favor of vectorized R, so the package no longer needs compilation and `bw_dir_lscv()` drops its `R_code` argument. With the default kernel, `bw_dir_lscv()` is 6--10 times faster. Both `bw_dir_lscv()` and `bw_dir_lcv()` keep memory flat in `n`, and `bw_dir_lcv()` no longer allocates a large three-dimensional array for non-von Mises--Fisher kernels.
* Fix `bw_dir_lscv()` ignoring `L`, which gave the von Mises--Fisher loss for any kernel.
* Fix `d_mixvmf()`, `R_Psi_mixvmf()`, `kappa_ml()`, and `bw_dir_rot()` for large concentrations. `kappa_ml()` now defaults to `max_kappa = 1e4`, and `bw_dir_rot()`, `bw_dir_ami()`, and `bw_dir_emi()` warn and return `NA` for degenerate concentrations.
* Add `sim` and `f_sim` to `int_hypsph()` for importance sampling, which fixes `bw_dir_emi()`, `bw_dir_ami()`, and `bw_dir_lscv()` for concentrated data with `q > 2`. `bw_dir_emi()` also computes its MISE integrand in log scale, avoiding over/underflow and loss of precision for small bandwidths. Document its overflow for large `q`.
* Optimize `bw_dir_lcv()`, `bw_dir_lscv()`, and `bw_dir_emi()` on `log(h)`, warning if `optim()` does not converge or stops at `optim_lower` or `optim_upper`. `optim_lower` now defaults to `1e-4` in `bw_dir_lcv()` and to `max(1e-4, exp(-700 / q) / sqrt(2 * pi))` in `bw_dir_lscv()` and `bw_dir_emi()`, keeping their losses finite for large `q`.
* Make `R_Psi_mixvmf()` stop instead of returning `NA` when its integration fails.
* Add `kappa_max` to `bic_vmf_mix()` to discard mixtures with a collapsed component. When no mixture can be fitted, `bic_vmf_mix()` now warns and returns the single von Mises--Fisher fit; a refined fit that worsens `crit` is discarded.
* Fix `to_rad()` mapping the S^2 south pole to the wrong polar angle, and `int_hypsph()` ignoring `na.rm`.
* Fix documentation and add `\link`/`\seealso` cross-references and a `testthat` unit-test suite.

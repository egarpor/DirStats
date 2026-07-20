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

# DirStats 0.2.0

* Fixed `to_rad()` mapping the \eqn{S^2} south pole to the wrong polar angle
  (`acos()` result was wrapped with `%% pi`, sending \eqn{\pi} to \eqn{0}).
* `int_hypsph()` now honors its `na.rm` argument (it was hardcoded to `TRUE`).
* `bic_vmf_mix()` now supports `crit = "AICc"` (previously errored).
* Efficiency: `kde_dir()` no longer computes the normalizing constant twice
  for user-defined kernels; `bw_dir_lcv()` avoids allocating a large
  three-dimensional array in the non-von Mises--Fisher branch.
* Documentation fixes and added `\link`/`\seealso` cross-references.
* Added a `testthat` unit-test suite.

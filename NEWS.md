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

* Fixed `to_rad()` mapping the \eqn{S^2} south pole to the wrong polar angle
  (`acos()` result was wrapped with `%% pi`, sending \eqn{\pi} to \eqn{0}).
* `int_hypsph()` now honors its `na.rm` argument (it was hardcoded to `TRUE`).
* `bic_vmf_mix()` no longer errors with "argument is of length zero" when no
  mixture can be fitted: the iterative search is skipped and the single von
  Mises--Fisher fallback is reached, now warning and returning `mu_hat` as a
  one-row matrix so that the plug-in selectors can consume it. The refined fit
  of the optimal mixture is also discarded if it worsens `crit`.
* `bic_vmf_mix()` gains `kappa_max` (defaults to `1e4`), which discards mixtures
  holding a component above that concentration. The likelihood of a mixture of
  von Mises--Fisher densities is unbounded as a component collapses onto a few
  observations, and `crit` then favors it: such a spike drove `bw_dir_ami()` an
  order of magnitude below `bw_dir_rot()`. Set `kappa_max = Inf` to disable the
  screening.
* Numerical robustness: the `q = 1` branch of `bw_dir_rot()` now uses
  exponentially scaled Bessel functions (consistent with the `q > 2` branch),
  giving identical values while avoiding overflow for large concentrations.
  The MISE integrand of `bw_dir_emi()` builds the expectation term in log
  scale, so it no longer over/underflows for small bandwidths (previously the
  `NaN`s were silently masked by the integration routines' `na.rm`).
* Efficiency: `kde_dir()` no longer computes the normalizing constant twice
  for user-defined kernels; `bw_dir_lcv()` avoids allocating a large
  three-dimensional array in the non-von Mises--Fisher branch.
* Removed the FORTRAN subroutines `lcv_dir_vmf` and `lscv_dir_vmf`, together
  with the vendored Netlib `ribesl` they relied on, in favor of vectorized R
  expressions built on `c_h()` and `log_c_vmf()`. The least-squares
  cross-validation is 9--10 times faster (78.8 to 8.3 seconds for `n = 2000`
  over 20 bandwidths), because `ribesl` had to be called once per pair and
  bandwidth, and the likelihood cross-validation is 1.6 times faster on the
  default optimization path, as the bandwidth-independent terms are now
  computed once. The pairwise sums of both selectors are accumulated in chunks
  of at most 2.5e5 inner products, so the memory footprint stays flat in `n`
  (`bw_dir_lcv()` caches the whole `n x n` matrix for `n <= 5000`, where
  recomputing it per evaluation costs more than the exponentials it feeds).
  Losses are unchanged to 3e-11, and were checked against the independent
  implementations in `polykde::bw_cv_polysph()`.
* `int_hypsph()` gains `sim` and `f_sim`, which turn it into an importance
  sampling integrator. Uniform draws spread the Monte Carlo effort evenly over
  the sphere and therefore miss integrands whose mass concentrates on a small
  region, collapsing the estimate to (almost) zero; the effect worsens with the
  dimension. `sim` accepts either a sampler function or a matrix of draws, the
  latter reusing the same sample across integrations. The default remains a
  uniform sample, with results unchanged bit for bit.
* Fixed `bw_dir_emi()`, `bw_dir_ami()`, and `bw_dir_lscv()` returning
  meaningless bandwidths on concentrated data for `q > 2`. All three integrate
  functionals of a mixture of von Mises--Fisher densities over the sphere, which
  the uniform Monte Carlo integration missed: the squared bias of `bw_dir_emi()`
  silently vanished from the MISE, leaving an objective with no interior
  minimum, and the curvature of `bw_dir_ami()` collapsed to zero. They now
  sample half the points from the reference mixture and half uniformly, the
  uniform part bounding the importance weights where the mixture is negligible
  but the integrand is not.
* Fixed `bw_dir_emi()`, `bw_dir_ami()`, and `bw_dir_lscv()` hanging forever
  (observed for more than 10 hours on a single call) when the reference mixture
  holds a collapsed component. The importance sampling proposal is drawn with
  `movMF::rmovMF()`, whose rejection sampler has no iteration ceiling and never
  returns for concentrations of the order of `1 / .Machine$double.eps`. The
  concentrations used to build the proposal are now capped at `1e10`, with a
  warning. Only the proposal is capped, which merely widens it and leaves the
  importance sampling estimator unbiased; the integrand keeps the concentrations
  of the mixture. Note that an infinite concentration errors instead of hanging,
  so it is the large-but-finite range that has to be bounded.
* `bw_dir_emi()` draws its Monte Carlo sample once instead of once per
  evaluation of the objective, so the MISE is a deterministic function of the
  bandwidth and `optim()` no longer differences noise when estimating its
  gradient. On a concentrated mixture on the 10-sphere, the selected bandwidth
  is now stable to about 1% across seeds, against a fivefold spread before.
* Fixed `bw_dir_lscv()` ignoring `L`: it was not passed down to the internal
  `kde_dir()` calls, so any kernel silently gave the von Mises--Fisher loss.
  Its leave-one-out terms are also vectorized now, over 10 times faster.
* Removed the `R_code` argument of `bw_dir_lscv()`, as there is no longer a
  FORTRAN path to contrast with the R one.
* New exported utilities `log_besselI_scaled()` and `log_c_vmf()` for the
  stable evaluation of the scaled log-Bessel function and the (log) von
  Mises--Fisher normalizing constant. They combine a cubic-spline table (orders
  up to 50), a small-argument series, and large-argument / large-order
  (DLMF 10.41) asymptotics, staying finite and accurate where `besselI()` would
  over/underflow or lose precision. Adds a dependency on the `Bessel` package.
* `kappa_ml()` now inverts the Bessel-ratio likelihood equation over
  `log(kappa)` using `log_besselI_scaled()`, with an auto-extending search
  bracket. This removes the previous hard `[1e-4, 100]` cap (highly concentrated
  data no longer errors), is stable in high dimensions, and matches
  `movMF::movMF()` to several significant figures. `max_kappa` now defaults to
  `1e4` and, together with `min_kappa`, is only an initial search range.
* `bw_dir_rot()` now evaluates the rule of thumb of Proposition 2 in
  García-Portugués (2013) with a single log-scale formula for all `q`, built on
  `log_besselI_scaled()` and a log-sum-exp, replacing the separate
  `q == 1` / `q == 2` / `q > 2` closed forms. This fixes the selector returning
  `NaN` for very concentrated data (concentrations above about `1e5`, now
  reachable since `kappa_ml()` is uncapped), because `besselI()` silently
  returns `0` there. Values are unchanged (to ~1e-11) for moderate
  concentrations. Degenerate concentration estimates now raise a warning and
  return `NA` instead of a silent `NaN`.
* `c_h()` is computed in log scale for the von Mises--Fisher kernel. It
  previously returned `Inf` for small bandwidths (roughly `h < 0.002`), because
  `besselI(1 / h^2, ...)` silently returns `0`; `kde_dir()` then failed with
  "NA/NaN/Inf in foreign function call". Such bandwidths are returned by
  `bw_dir_rot()` for concentrated data, so the selector and the estimator are
  now coherent.
* The exact expression of `R_Psi_mixvmf()` is computed in log scale, fixing the
  curvature collapsing to `0` (concentrations around `1e5`) or `NaN` (around
  `1e6`). A null curvature was particularly harmful, as it appears in a
  denominator of the `bw_dir_ami()` bandwidth.
* `d_mixvmf()` and the numeric branch of `R_Psi_mixvmf()` use
  `log_c_vmf()` instead of `rotasym::c_vMF()`, which returns `Inf` for
  concentrations above about `1e5`. Mixtures holding a concentrated component
  previously made `bw_dir_emi()` return `Inf` and `bw_dir_ami()` return `0`.
  Values are unchanged (to ~1e-10) for moderate concentrations.
* `log_besselI_scaled()` returns a zero-length vector for zero-length input,
  as `besselI()` does.
* Move `rotasym` from `Imports` to `Suggests`. The package no longer calls it at
  run time: the surface area of the hypersphere and the uniform sampling in
  `int_hypsph()` are now computed internally.
* Documentation fixes and added `\link`/`\seealso` cross-references.
* Added a `testthat` unit-test suite.

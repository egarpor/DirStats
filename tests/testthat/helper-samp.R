
# Draws on S^q obtained by projecting a N(sqrt(kappa) * mu, I_{q + 1}) onto the
# sphere, so that tests only needing spherical data with a given spread do not
# depend on the suggested rotasym. Scaling the mean by sqrt(kappa) matches the
# vMF spread, as 1 - E[X'mu] is q / (2 * kappa) for both
r_proj_norm <- function(n, mu, kappa) {

  p <- length(mu)
  return(normalize(matrix(rnorm(n = n * p, mean = sqrt(kappa) * mu), nrow = n,
                          ncol = p, byrow = TRUE)))

}

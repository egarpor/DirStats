

#' @title \code{DirStats} -- Nonparametric Methods for Directional Data
#'
#' @description Nonparametric kernel density estimation, bandwidth selection,
#' and other utilities for analyzing directional data. Implements the estimator
#' in Bai, Rao and Zhao (1987) <doi:10.1016/0047-259X(88)90113-3>, the
#' cross-validation bandwidth selectors in Hall, Watson and Cabrera (1987)
#' <doi:10.1093/biomet/74.4.751> and the plug-in bandwidth selectors in
#' García-Portugués (2013) <doi:10.1214/13-ejs821>.
#'
#' @author Eduardo García-Portugués.
#' @docType package
#' @name DirStats-package
#' @import graphics stats utils grDevices
#' @useDynLib DirStats, .registration = TRUE
#' @aliases DirStats DirStats-package
NULL

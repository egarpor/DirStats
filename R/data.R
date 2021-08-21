

#' @title Lebedev quadrature on the sphere
#'
#' @description Nodes and weights for Lebedev quadrature on the sphere
#' \eqn{S^2}. The rule has 5810 points and is exact up to polynomials
#' of order 131.
#'
#' @docType data
#' @format A data frame with 5810 rows and two variables:
#' \describe{
#'   \item{xyz}{nodes for quadrature, a matrix with three columns.}
#'   \item{w}{weights for quadrature, a vector.}
#' }
#' @details
#' The approximation to the integral of \eqn{f} has the form
#' \deqn{\int_{S^2} f(x, y, z) \,\mathrm{d}x \,\mathrm{d}y \,\mathrm{d}z =
#' 4 \pi \sum_{i = 1}^N w_i f(x_i, y_i, z_i)}{\int_{S^2} f(x, y, z)
#' dx dy dz = 4 \pi \sum_{i = 1}^N w_i f(x_i, y_i, z_i)}
#' where \eqn{N = 5810}. The nodes (in spherical coordinates) and weights
#' are processed from
#' \href{https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/lebedev_131.txt
#' }{lebedev_131.txt}.
#' @source \url{
#' https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html
#' }
#' @references
#' Lebedev, V. I. and Laikov, D. N. (1999). A quadrature formula for the
#' sphere of the 131st algebraic order of accuracy. \emph{Doklady Mathematics},
#' 59(3):477--481.
#' @examples
#' # Load data
#' data("lebedev")
#'
#' # Integrate x_1 * x_2^2 (zero integral)
#' f_1 <- function(x) x[, 1] * x[, 2]^2
#' 4 * pi * sum(lebedev$w * f_1(lebedev$xyz))
"lebedev"

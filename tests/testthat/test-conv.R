
## norm2()

test_that("norm2 works for vectors and matrices", {

  expect_equal(norm2(c(3, 4)), 5)
  expect_equal(norm2(rbind(c(3, 4), c(0, 1), c(1, 0))), c(5, 1, 1))

})

## normalize()

test_that("normalize returns unit-norm rows", {

  x <- rbind(1:3, c(-2, 5, 1), c(10, 0, 0))
  expect_equal(norm2(normalize(x)), rep(1, nrow(x)))
  expect_equal(norm2(normalize(1:4)), 1)

})

## to_cir() and to_rad()

test_that("to_cir and to_rad are inverses on the circle", {

  th <- seq(0, 2 * pi, l = 20)[-20]
  expect_equal(to_rad(to_cir(th)), th)
  expect_equal(drop(to_cir(0)), c(1, 0))

})

test_that("to_rad maps the poles of S^2 correctly (regression: south pole)", {

  expect_equal(to_rad(c(0, 0, 1))[2], 0)
  expect_equal(to_rad(c(0, 0, -1))[2], pi)
  r <- to_rad(c(0, 0, -1))
  expect_equal(drop(to_sph(r[1], r[2])), c(0, 0, -1))

})

test_that("to_rad only accepts circular or spherical data", {

  expect_error(to_rad(diag(4)), "circular or spherical")

})

## to_sph() and to_rad()

test_that("to_sph and to_rad are inverses on the sphere", {

  th <- c(0.3, 1.2, 4.5)
  ph <- c(0.1, 1.5, 2.9)
  rad <- to_rad(to_sph(th, ph))
  expect_equal(rad[, 1], th)
  expect_equal(rad[, 2], ph)

})

## w_q()

test_that("w_q matches rotasym::w_p and its log", {

  skip_if_not_installed("rotasym")
  q <- 0:10
  expect_equal(w_q(q = q), rotasym::w_p(p = q + 1))
  expect_equal(w_q(q = q, log = TRUE), log(w_q(q = q)))

})

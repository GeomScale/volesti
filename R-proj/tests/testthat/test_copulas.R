context("Copulas' test")

library(volesti)

test_that("10-dimensional 2-hyp_fam copula", {
  h1 = h1 = runif(n = 10, min = 1, max = 1000)
  h1 = h1 / 1000
  h2 = h2 = runif(n = 10, min = 1, max = 1000)
  h2 = h1 / 1000
  cop = copula1(h1 = h1 , h2 = h2 , numSlices = 10, N = 100000)
  res = sum(cop)
  expect_equal(res, 1)
})

test_that("20-dimensional 1_hyp_1_ell fam copula", {
  h = h = runif(n = 20, min = 1, max = 1000)
  h = h / 1000
  E = replicate(20, rnorm(30))
  E = cov(E)
  cop = copula2(h = h , E = E , numSlices = 100, N = 100000)
  res = sum(cop)
  expect_equal(res, 1)
})
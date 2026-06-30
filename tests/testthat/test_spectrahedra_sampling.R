context("Spectrahedron sampling")

# These tests validate the sampling functionality for Spectrahedron objects
# using various random walk types. The tests check basic functionality,
# parameter validation, and sampling quality metrics.

test_that("sample_points method exists for Spectrahedron", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  
  expect_true(isGeneric("sample_points"))
  expect_true(hasMethod("sample_points", "Spectrahedron"))
})

test_that("sample_points requires valid Spectrahedron object", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  
  expect_error(sample_points(list(), N = 100), "Spectrahedron only")
})

test_that("sample_points validates parameters", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  # Invalid number of samples (N)
  expect_error(sample_points(spec, N = 0), 
               "N must be a positive integer")
  
  expect_error(sample_points(spec, N = -1), 
               "N must be a positive integer")
  
  # Invalid walk type
  expect_error(sample_points(spec, N = 100, walk_type = "INVALID"), 
               "walk_type must be one of")
  
  # Invalid number of burns
  expect_error(sample_points(spec, N = 100, n_burns = -1), 
               "n_burns must be non-negative")
})

test_that("sample_points returns correct dimensions", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  d <- dim(spec)
  N <- 50
  
  samples <- sample_points(spec, N = N, walk_type = "RDHR")
  
  expect_is(samples, "matrix")
  expect_equal(nrow(samples), d)
  expect_equal(ncol(samples), N)
})

test_that("sample_points works with all supported walk types", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  walk_types <- c("RDHR", "CDHR", "HMC", "BILLIARD")
  
  for (walk in walk_types) {
    samples <- sample_points(spec, N = 30, walk_type = walk)
    
    expect_is(samples, "matrix")
    expect_equal(nrow(samples), dim(spec))
    expect_equal(ncol(samples), 30)
  }
})

test_that("sample_points respects random seed", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  set.seed(42)
  samples1 <- sample_points(spec, N = 40, walk_type = "RDHR")
  
  set.seed(42)
  samples2 <- sample_points(spec, N = 40, walk_type = "RDHR")
  
  # Should be identical with same seed
  expect_equal(samples1, samples2)
})

test_that("sample_points works with burn-in steps", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  set.seed(123)
  samples_no_burn <- sample_points(spec, N = 50, walk_type = "RDHR", 
                                   n_burns = 0)
  
  set.seed(123)
  samples_burn <- sample_points(spec, N = 50, walk_type = "RDHR", 
                                n_burns = 100)
  
  # Both should have same dimensions
  expect_equal(dim(samples_no_burn), dim(samples_burn))
  
  # Samples should differ due to burn-in
  expect_false(identical(samples_no_burn, samples_burn))
})

test_that("sample_points generates samples in the feasible region", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  samples <- sample_points(spec, N = 100, walk_type = "RDHR")
  
  # Check that each sample is feasible
  # (All samples should pass the feasibility check)
  for (i in 1:ncol(samples)) {
    point <- samples[, i]
    in_spec <- is_point_in_spectrahedron(spec, point)
    expect_true(in_spec, info = paste("Sample", i, "is not feasible"))
  }
})

test_that("sample_points with short burn-in", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  samples <- sample_points(spec, N = 30, walk_type = "RDHR", n_burns = 10)
  
  expect_is(samples, "matrix")
  expect_equal(nrow(samples), dim(spec))
  expect_equal(ncol(samples), 30)
})

test_that("sample_points with HMC walk", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  samples <- sample_points(spec, N = 20, walk_type = "HMC")
  
  expect_is(samples, "matrix")
  expect_equal(ncol(samples), 20)
})

test_that("sample_points works with default parameters", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  # Should work with only N parameter
  samples <- sample_points(spec, N = 50)
  
  expect_is(samples, "matrix")
  expect_equal(nrow(samples), dim(spec))
  expect_equal(ncol(samples), 50)
})

test_that("sample_points produces sufficient samples for analysis", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  N <- 200
  samples <- sample_points(spec, N = N, walk_type = "RDHR")
  
  # Verify sample size and dimensions are as expected
  expect_equal(ncol(samples), N)
  expect_equal(nrow(samples), dim(spec))
  
  # Check numerical properties
  expect_true(all(is.finite(samples)))
  expect_false(any(is.na(samples)))
})

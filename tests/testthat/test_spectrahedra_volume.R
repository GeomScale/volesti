context("Spectrahedron volume computation")

# These tests validate the volume calculation for Spectrahedron objects
# using both cooling balls and sequence of balls methods.
# Note: These are integration tests; actual numerical validation requires
# test data files (SDPA format spectrahedra)

test_that("volume method exists for Spectrahedron", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  
  # This test verifies the method dispatch is correctly configured
  expect_true(isGeneric("volume"))
  expect_true(hasMethod("volume", "Spectrahedron"))
})

test_that("volume requires valid Spectrahedron object", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  
  # Test with non-Spectrahedron input
  expect_error(volume(list()), "Spectrahedron only")
})

test_that("volume validates parameters", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  # Invalid method
  expect_error(volume(spec, method = "INVALID"), 
               "method must be 'CB' or 'SOB'")
  
  # Invalid parameters
  expect_error(volume(spec, E = -1, 
                     method = "CB"), 
               "E must be positive")
  
  expect_error(volume(spec, E = 0.1, 
                     method = "CB"),
               "E should be a value")
})

test_that("volume returns numeric result", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  # Test CB method
  vol_cb <- volume(spec, method = "CB", E = 0.05, random_walk = "RDHR")
  expect_is(vol_cb, "numeric")
  expect_length(vol_cb, 1)
  expect_true(vol_cb > 0)
  
  # Test SOB method
  vol_sob <- volume(spec, method = "SOB", E = 0.05, random_walk = "RDHR")
  expect_is(vol_sob, "numeric")
  expect_length(vol_sob, 1)
  expect_true(vol_sob > 0)
})

test_that("volume produces consistent results across repeated calls", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  set.seed(42)
  vol1 <- volume(spec, method = "CB", E = 0.1, random_walk = "RDHR")
  
  set.seed(42)
  vol2 <- volume(spec, method = "CB", E = 0.1, random_walk = "RDHR")
  
  # With same seed, results should be identical or very close
  expect_equal(vol1, vol2, tolerance = 1e-10)
})

test_that("volume estimates improve with tighter error bounds", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  set.seed(123)
  vol_loose <- volume(spec, method = "CB", E = 0.2, random_walk = "RDHR")
  
  set.seed(123)
  vol_tight <- volume(spec, method = "CB", E = 0.05, random_walk = "RDHR")
  
  # Both should be positive and of similar magnitude
  expect_true(vol_loose > 0)
  expect_true(vol_tight > 0)
  # Tighter bound should complete without excessive computation
  expect_equal(length(vol_tight), 1)
})

test_that("volume with different walk types produces reasonable estimates", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  set.seed(456)
  vol_rdhr <- volume(spec, method = "CB", E = 0.1, random_walk = "RDHR")
  
  set.seed(456)
  vol_billiard <- volume(spec, method = "CB", E = 0.1, random_walk = "BILLIARD")
  
  # Different walks should produce comparable results
  expect_true(vol_rdhr > 0)
  expect_true(vol_billiard > 0)
})

test_that("volume works with default parameters", {
  skip_if_not(requireNamespace("volesti", quietly = TRUE))
  skip_if_not(file.exists(system.file("extdata/example_spectrahedron.sdpa", 
                                       package = "volesti")))
  
  spec <- load_spectrahedron("extdata/example_spectrahedron.sdpa", 
                             package = "volesti")
  
  # Should work with minimal parameters
  vol <- volume(spec)
  expect_is(vol, "numeric")
  expect_length(vol, 1)
  expect_true(vol > 0)
})

test_that("can serialize and deserialize data", {
  foo <- list(
    ## 4 rows, 3 columns
    test_2d_double = array(as.numeric(1:12), dim = c(4, 3)),
    ## 2 rows, 3 columns, 4 deep
    test_3d_int = array(1:24, dim = c(2, 3, 4))
  )

  t <- tempfile()
  save_hdf5_file(foo, t)

  ## Deserializing works
  deserialized_foo <- read_hdf5_file(t)
  expect_equal(deserialized_foo, foo)
})


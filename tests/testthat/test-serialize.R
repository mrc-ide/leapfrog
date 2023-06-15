test_that("can serialize and deserialize eigen data", {
  foo <- list(
    ## 4 rows, 3 columns
    test_2d_double = array(as.numeric(1:12), dim = c(4, 3)),
    ## 2 rows, 3 columns, 4 deep
    test_3d_int = array(1:24, dim = c(2, 3, 4))
  )

  t1 <- tempfile()
  t2 <- tempfile()
  out <- serialize_vector(foo, t1, t2)

  expect_equal(out, list(
    test_2d_double_path = t1,
    test_3d_int_path = t2
  ))

  content <- readLines(out$test_2d_double_path)
  expect_length(content, 3)
  expect_equal(content[1], "double")
  expect_equal(content[2], "4,3")
  expect_equal(content[3], paste0(1:12, collapse = ","))

  content <- readLines(out$test_3d_int_path)
    expect_length(content, 3)
    expect_equal(content[1], "int")
    expect_equal(content[2], "2,3,4")
    expect_equal(content[3], paste0(1:24, collapse = ","))

  ## Deserializing works
  deserialized <- deserialize_vector(t1, t2)
  expect_equal(deserialized, foo)
})

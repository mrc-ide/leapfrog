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

test_that("can serialize from R", {
  foo <- list(
    ## 4 rows, 3 columns
    test_2d_double = array(as.numeric(1:12), dim = c(4, 3)),
    ## 2 rows, 3 columns, 4 deep
    test_3d_int = array(1:24, dim = c(2, 3, 4))
  )

  t1 <- tempfile()
  t2 <- tempfile()
  out1 <- serialize_r_to_tensor(foo$test_2d_double, t1)
  out2 <- serialize_r_to_tensor(foo$test_3d_int, t2)

  content <- readLines(out1)
  expect_length(content, 3)
  expect_equal(content[1], "double")
  expect_equal(content[2], "4,3")
  expect_equal(content[3], paste0(1:12, collapse = ","))

  content <- readLines(out2)
  expect_length(content, 3)
  expect_equal(content[1], "int")
  expect_equal(content[2], "2,3,4")
  expect_equal(content[3], paste0(1:24, collapse = ","))

  ## Deserializing works
  deserialized <- deserialize_vector(out1, out2)
  expect_equal(deserialized, foo)
})

test_that("can deseralize from R", {
  foo <- list(
    ## 4 rows, 3 columns
    test_2d_double = array(as.numeric(1:12), dim = c(4, 3)),
    ## 2 rows, 3 columns, 4 deep
    test_3d_int = array(1:24, dim = c(2, 3, 4))
  )

  t1 <- tempfile()
  t2 <- tempfile()
  out1 <- serialize_r_to_tensor(foo$test_2d_double, t1)
  out2 <- serialize_r_to_tensor(foo$test_3d_int, t2)

  d1 <- deserialize_tensor_to_r(out1)
  d2 <- deserialize_tensor_to_r(out2)

  expect_identical(foo$test_2d_double, d1)
  expect_identical(foo$test_3d_int, d2)
})

test_that("R can serialize 1d vector", {
  data <- array(c(1L, 2L, 3L, 10L), dim = 4)
  t <- tempfile()
  path <- serialize_r_to_tensor(data, t)
  content <- readLines(path)
  expect_length(content, 3)
  expect_equal(content[1], "int")
  expect_equal(content[2], "4")
  expect_equal(content[3], paste0(data, collapse = ","))

  d <- deserialize_tensor_to_r(path)

  expect_equal(d, data)
})

test_that("R can serialize boolean type as int", {
  data <- array(c(TRUE, FALSE, FALSE, TRUE), dim = c(2, 2))
  t <- tempfile()
  path <- serialize_r_to_tensor(data, t)
  content <- readLines(path)
  expect_length(content, 3)
  expect_equal(content[1], "int")
  expect_equal(content[2], "2,2")
  expect_equal(content[3], paste0(c(1, 0, 0, 1), collapse = ","))
})

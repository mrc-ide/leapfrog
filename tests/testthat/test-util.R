test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})


test_that("can format a vector", {
  expect_equal(format_vector("thing"), "'thing'")
  expect_equal(format_vector(c("foo", "bar")), "'foo', 'bar'")
})

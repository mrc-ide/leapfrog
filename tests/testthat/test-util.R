test_that("null-or-value works", {
  expect_equal(1 %||% NULL, 1)
  expect_equal(1 %||% 2, 1)
  expect_equal(NULL %||% NULL, NULL)
  expect_equal(NULL %||% 2, 2)
})

test_that("can assert value is from valid set", {
  expect_true(assert_enum("foo", c("foo", "bar")))
  input <- "thing"
  expect_error(assert_enum(input, c("foo", "bar")),
               "input must be one of 'foo', 'bar', got 'thing'")
})

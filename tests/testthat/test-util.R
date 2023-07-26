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

test_that("assert_single_item_set", {
  input <- list(foo = "1", bar = "2")
  expect_silent(assert_single_item_set(input, "foo"))
  expect_silent(assert_single_item_set(input, "bar"))
  expect_error(
    assert_single_item_set(input, c("foo", "bar")),
    "Items 'foo', 'bar' all have values, only one should be set",
    fixed = TRUE)
})


test_that("can format a vector", {
  expect_equal(format_vector("thing"), "'thing'")
  expect_equal(format_vector(c("foo", "bar")), c("'foo', 'bar'"))
})

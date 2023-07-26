test_that("assert_scalar", {
  expect_error(assert_scalar(NULL), "must be a scalar")
  expect_error(assert_scalar(numeric(0)), "must be a scalar")
  expect_error(assert_scalar(1:2), "must be a scalar")
})


test_that("assert_character", {
  expect_silent(assert_character("a"))
  expect_error(assert_character(1), "must be character")
  expect_error(assert_character(TRUE), "must be character")
})

test_that("can assert value is from valid set", {
  expect_silent(assert_enum("foo", c("foo", "bar")))
  input <- "thing"
  expect_error(assert_enum(input, c("foo", "bar")),
               "input must be one of 'foo', 'bar', got 'thing'")
})

test_that("assert_one_is_set", {
  input <- list(foo = "1", bar = "2")
  expect_silent(assert_one_is_set(input, "foo"))
  expect_silent(assert_one_is_set(input, "bar"))
  expect_error(
    assert_one_is_set(input, c("foo", "bar")),
    "Items 'foo', 'bar' all have values, only one should be set",
    fixed = TRUE)
})

test_that("assert_set", {
  expect_silent(assert_set("value"))
  expect_error(assert_set("",
                          "\"\" must be a non-null, non-na, non-empty string"))
  expect_error(assert_set(NA_character_,
                          "\"\" must be a non-null, non-na, non-empty string"))
  expect_error(assert_set(NULL,
                          "\"\" must be a non-null, non-na, non-empty string"))
  expect_error(assert_set("  ",
                          "\"\" must be a non-null, non-na, non-empty string"))
})

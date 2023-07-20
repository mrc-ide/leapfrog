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

test_that("assert_names", {
  required <- c("one", "two")
  optional <- c("three", "four")
  input <- list(one = 1, two = 2, three = 3, four = 4)
  expect_true(assert_names(input, required, optional))

  input <- list(one = 1, two = 2)
  expect_true(assert_names(input, required, optional))

  input <- list(one = 1, three = 3, four = 4)
  expect_error(assert_names(input, required, optional),
               "Required item(s) 'two' are missing from input",
               fixed = TRUE)

  input <- list(one = 1, two = 2, five = 5)
  expect_error(assert_names(input, required, optional),
               "Unknown item(s) 'five' are included in input",
               fixed = TRUE)
})

test_that("assert_names_one_of", {
  input <- list(foo = "1", bar = "2")
  assert_names_one_of(input, c("foo"))
  assert_names_one_of(input, c("bar"))
  expect_error(
    assert_names_one_of(input, c("foo", "bar")),
    "Items 'foo', 'bar' are all included in input, only one should be present",
    fixed = TRUE)
})

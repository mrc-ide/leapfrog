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


test_that("can group list of lists", {
  input_list <- list(
    list(id = 1, type = "foo"),
    list(id = 2, type = "bar"),
    list(id = 3, type = "foo"),
    list(id = 4, type = "bar")
  )
  expect_equal(group_list_of_lists(input_list, "type"),
               list(
                 "foo" = list(
                   list(id = 1, type = "foo"),
                   list(id = 3, type = "foo")
                 ),
                 "bar" = list(
                   list(id = 2, type = "bar"),
                   list(id = 4, type = "bar")
                 )
               ))
})

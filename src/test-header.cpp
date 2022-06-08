#include <testthat.h>

#include "code.h"

context("C++ code can be unit tested") {

  test_that("header works") {
    expect_true(add(2, 2) == 4);
  }
}

#include <Rcpp.h>

#include "code.hpp"

//' Add two numbers.
//'
//' @param a First number.
//' @param b Second number.
//' @return The sum of two values.
//' @export
// [[Rcpp::export]]
int adder(int a, int b) {
  return add(a, b);
}

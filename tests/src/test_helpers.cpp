#include <Rcpp.h>

#include "frogger.hpp"
#include "types.hpp"
#include "serialize_eigen.hpp"

// [[Rcpp::export]]
Rcpp::List serialize_vector(const Rcpp::List data, const std::string path1, const std::string path2) {
  // 4 rows, 3 columns
  const leapfrog::TensorMap2<double> test_2d(REAL(data["test_2d_double"]), 4, 3);
  const leapfrog::TensorMap3<int> test_3d(INTEGER(data["test_3d_int"]), 2, 3, 4);

  serialize::serialize_tensor_map(test_2d, path1);
  serialize::serialize_tensor_map(test_3d, path2);

  Rcpp::List ret =
      Rcpp::List::create(Rcpp::_["test_2d_double_path"] = path1,
                         Rcpp::_["test_3d_int_path"] = path2);
  return ret;
}

// [[Rcpp::export]]
Rcpp::List deserialize_vector(const std::string path1, const std::string path2) {
  const Eigen::Tensor<double, 2> test_2d = serialize::deserialize_tensor<double, 2>(path1);
  const Eigen::Tensor<int, 3> test_3d = serialize::deserialize_tensor<int, 3>(path2);

  Rcpp::NumericVector r_test_2d(4 * 3);
  Rcpp::IntegerVector r_test_3d(2 * 3 * 4);

  r_test_2d.attr("dim") =
      Rcpp::NumericVector::create(4, 3);
  r_test_3d.attr("dim") =
      Rcpp::NumericVector::create(2, 3, 4);

  std::copy_n(test_2d.data(), test_2d.size(), REAL(r_test_2d));
  std::copy_n(test_3d.data(), test_3d.size(), INTEGER(r_test_3d));

  Rcpp::List ret =
      Rcpp::List::create(Rcpp::_["test_2d_double"] = r_test_2d,
                         Rcpp::_["test_3d_int"] = r_test_3d);
  return ret;
}

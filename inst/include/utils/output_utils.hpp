#pragma once

#include "language_types.hpp"
#include "../generated/concepts.hpp"

#ifdef CALLER_R
#include <Rcpp.h>

namespace leapfrog {
namespace internal {

template<typename TensorType, typename... Args>
void write_output_data(OutputData &ret, const TensorType &tensor, const std::string name, int idx, Args... dims) {
  int size = (dims * ... * 1);
  Rcpp::NumericVector r_vec(size);
  r_vec.attr("dim") = Rcpp::IntegerVector::create(dims...);
  std::copy_n(tensor.data(), tensor.size(), REAL(r_vec));
  ret.names[idx] = name;
  ret.data[idx] = r_vec;
}

}
}
#endif

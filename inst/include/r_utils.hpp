#pragma once

#include <Rcpp.h>
#include <unsupported/Eigen/CXX11/Tensor>

template <typename T>
T* r_data(SEXP x) {
  static_assert(sizeof(T) == 0, "Only specializations of r_data can be used");
}

template <>
double* r_data(SEXP x) {
  return REAL(x);
}

template <>
int * r_data(SEXP x) {
  return INTEGER(x);
}

template<typename T, std::size_t rank>
auto convert_0_based(const Eigen::TensorMap<Eigen::Tensor<int, rank>> map) {
  static_assert(sizeof(T) == 0, "Only specializations of convert_0_based can be used");
}

template<std::size_t rank>
auto convert_0_based(const Eigen::TensorMap<Eigen::Tensor<double, rank>> map) {
  Eigen::Tensor<double, rank> new_tensor = map; // Create a copy
  for (int i = 0; i < new_tensor.size(); ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    new_tensor.data()[i] = new_tensor.data()[i] - 1.0f;
  }
  return new_tensor;
}

template<std::size_t rank>
auto convert_0_based(const Eigen::TensorMap<Eigen::Tensor<int, rank>> map) {
  Eigen::Tensor<int, rank> new_tensor = map; // Create a copy
  for (int i = 0; i < new_tensor.size(); ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    new_tensor.data()[i] = new_tensor.data()[i] - 1;
  }
  return new_tensor;
}

auto convert_0_based(std::vector<int>& input) {
  for (auto& num: input) {
    num -= 1;
  }
}

template<typename T, typename... Args>
auto parse_data(const Rcpp::List data, const std::string& key, Args... dims) {
  constexpr std::size_t rank = sizeof...(dims);
  Eigen::array<int, rank> dimensions{ static_cast<int>(dims)... };

  int length = std::accumulate(dimensions.begin(), dimensions.end(), 1, std::multiplies<int>());
  SEXP array_data = data[key];
  // In cases where the input data has project years we might not use all of it model fit
  // So we can take create a Map over a smaller slice of the data
  // As long as this is true we can be confident we're not referencing invalid memory
  if (LENGTH(array_data) < length) {
    Rcpp::stop("Invalid size of data for '%s', expected %d got %d",
               key,
               length,
               LENGTH(array_data));
  }

  return Eigen::TensorMap<Eigen::Tensor<T, rank>>(r_data<T>(array_data), static_cast<int>(dims)...);
}

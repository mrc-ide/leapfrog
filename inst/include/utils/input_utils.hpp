#pragma once

#include "language_types.hpp"
#include <unsupported/Eigen/CXX11/Tensor>

#ifdef CALLER_R
#include <Rcpp.h>
#include "../generated/name_maps.hpp"

template <typename T>
T* r_data(SEXP x) {
  static_assert(sizeof(T) == 0, "Only specializations of r_data can be used");
}

template <>
double* r_data(SEXP x) {
  return REAL(x);
}

template <>
int* r_data(SEXP x) {
  return INTEGER(x);
}

template<typename T, typename... Args>
auto parse_data(const InputData &data, const std::string& cpp_name, Args... dims) {
  constexpr std::size_t rank = sizeof...(dims);
  // TODO: assert this name is valid and return good error if it is not
  const std::string& r_name = leapfrog::CPP_TO_R_NAME.at(cpp_name);
  if constexpr (rank == 0) {
      return Rcpp::as<T>(data[r_name]);
  } else {

    Eigen::array<int, rank> dimensions{ static_cast<int>(dims)... };

    int length = std::accumulate(dimensions.begin(), dimensions.end(), 1, std::multiplies<int>());
    SEXP array_data = data[r_name];
    // In cases where the input data has project years we might not use all of it model fit
    // So we can take create a Map over a smaller slice of the data
    // As long as this is true we can be confident we're not referencing invalid memory
    if (LENGTH(array_data) < length) {
        Rcpp::stop("Invalid size of data for '%s', expected %d got %d",
                r_name,
                length,
                LENGTH(array_data));
    }

    return Eigen::TensorMap<Eigen::Tensor<T, rank>>(r_data<T>(array_data), static_cast<int>(dims)...);
  }
}

#else

#include <filesystem>
#include "../serialize_eigen.hpp"

template<typename T, typename... Args>
auto parse_data(const InputData &input_dir, const std::string &key, Args... dims) {
  if (!std::filesystem::exists(input_dir / key)) {
    std::string msg = "File '" + key + "' in input dir '" + input_dir.string() + "' does not exist.\n";
    throw std::runtime_error(msg);
  }

  constexpr std::size_t rank = sizeof...(dims);
  if constexpr (rank == 0) {
    // TODO: impl
    return 1.0;
  } else {
    return serialize::deserialize_tensor<T, rank>(input_dir / key);
  }
}

#endif

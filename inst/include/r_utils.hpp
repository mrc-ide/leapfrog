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
auto convert_base(Eigen::TensorMap<Eigen::Tensor<int, rank>> map) {
  static_assert(sizeof(T) == 0, "Only specializations of convert_base can be used");
}

template<std::size_t rank>
auto convert_base(Eigen::TensorMap<Eigen::Tensor<double, rank>> map) {
  for (int i = 0; i < map.size(); ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    map.data()[i] = map.data()[i] - 1.0f;
  }
  return map;
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

template<std::size_t rank>
auto convert_base(Eigen::TensorMap<Eigen::Tensor<int, rank>> map) {
  for (int i = 0; i < map.size(); ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    map.data()[i] = map.data()[i] - 1;
  }
  return map;
}


// ListBuilder class copied from Kevin Ushey comment on Github issue:
// https://github.com/RcppCore/Rcpp/issues/243#issuecomment-73378636

// Usage example:
//
// // [[Rcpp::export]]
// List test_builder(SEXP x, SEXP y, SEXP z, std::string w) {
//   return ListBuilder()
//     .add("foo", x)
//     .add("bar", y)
//     .add("baz", z)
//     .add("bat", w);
// }
//
// /*** R
// test_builder(1:5, letters[1:5], rnorm(5), "abc")
// */

#include <Rcpp.h>
using namespace Rcpp;

class ListBuilder {

public:

  ListBuilder() {};
  ~ListBuilder() {};

  inline ListBuilder& add(const std::string& name, SEXP x) {
    names.push_back(name);
    elements.push_back(PROTECT(x));
    return *this;
  }

  template <typename T>
  inline ListBuilder& add(const std::string& name, const T& x) {
    names.push_back(name);
    elements.push_back(PROTECT(wrap(x)));
    return *this;
  }

  inline operator List() const {
    List result(elements.size());
    for (size_t i = 0; i < elements.size(); ++i) {
      result[i] = elements[i];
    }
    result.attr("names") = wrap(names);
    UNPROTECT(elements.size());
    return result;
  }

  inline operator DataFrame() const {
    List result = static_cast<List>(*this);
    result.attr("class") = "data.frame";
    result.attr("row.names") = IntegerVector::create(NA_INTEGER, XLENGTH(elements[0]));
    return result;
  }

private:

  std::vector<std::string> names;
  std::vector<SEXP> elements;

  ListBuilder(ListBuilder const&) {};

};

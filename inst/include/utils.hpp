#pragma once
#include <unsupported/Eigen/CXX11/Tensor>

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

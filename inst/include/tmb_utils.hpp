template<typename... Args>
auto convert_0_based_tmb_new(array<int> arr, Args... dims) {
  constexpr std::size_t rank = sizeof...(dims);
  Eigen::Tensor<int, rank> new_tensor(static_cast<int>(dims)...); // Create a copy

  for (int i = 0; i < arr.size(); ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    new_tensor.data()[i] = arr.data()[i] - 1;
  }
  return new_tensor;
}

template<typename Type, typename... Args>
auto convert_0_based_tmb_new(array<Type> arr, Args... dims) {
  constexpr std::size_t rank = sizeof...(dims);
  Eigen::Tensor<Type, rank> new_tensor(dims...); // Create a copy

  for (int i = 0; i < new_tensor.size(); ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    new_tensor.data()[i] = arr.data()[i] - 1.0f;
  }
  return new_tensor;
};

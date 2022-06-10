#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

template <typename Type>
using TensorMapX1cT = Eigen::TensorMap<Eigen::Tensor<const Type, 1>>;

template <typename Type>
using TensorMapX2cT = Eigen::TensorMap<Eigen::Tensor<const Type, 2>>;

template <typename Type>
using TensorMapX1T = Eigen::TensorMap<Eigen::Tensor<Type, 1>>;

template <typename Type>
using TensorMapX2T = Eigen::TensorMap<Eigen::Tensor<Type, 2>>;
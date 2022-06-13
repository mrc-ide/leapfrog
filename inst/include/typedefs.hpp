#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

template <typename real_type>
using TensorMapX1cT = Eigen::TensorMap<Eigen::Tensor<const real_type, 1>>;

template <typename real_type>
using TensorMapX2cT = Eigen::TensorMap<Eigen::Tensor<const real_type, 2>>;

template <typename real_type>
using TensorMapX1T = Eigen::TensorMap<Eigen::Tensor<real_type, 1>>;

template <typename real_type>
using TensorMapX2T = Eigen::TensorMap<Eigen::Tensor<real_type, 2>>;
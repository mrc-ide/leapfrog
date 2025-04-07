#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

namespace leapfrog {
namespace internal {

template<typename real_type>
using TM1 = Eigen::TensorMap<Eigen::Tensor<real_type, 1>>;

template<typename real_type>
using TM2 = Eigen::TensorMap<Eigen::Tensor<real_type, 2>>;

template<typename real_type>
using TM3 = Eigen::TensorMap<Eigen::Tensor<real_type, 3>>;

template<typename real_type>
using TM4 = Eigen::TensorMap<Eigen::Tensor<real_type, 4>>;


template<typename real_type>
using T1 = Eigen::Tensor<real_type, 1>;

template<typename real_type>
using T2 = Eigen::Tensor<real_type, 2>;

template<typename real_type>
using T3 = Eigen::Tensor<real_type, 3>;

template<typename real_type>
using T4 = Eigen::Tensor<real_type, 4>;

template<typename real_type>
using T5 = Eigen::Tensor<real_type, 5>;


template<typename real_type, typename std::ptrdiff_t... Dims>
using TFS = Eigen::TensorFixedSize<real_type, Eigen::Sizes<Dims...>>;

}
}

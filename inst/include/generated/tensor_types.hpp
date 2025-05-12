#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

namespace leapfrog {
namespace internal {

#ifdef CALLER_CPP

// When we run the model from the "CPP" interface we read input data from disk
// when we run it from R we map over data owned by R. So from CPP caller we need
// to operate on owned data, not the mapped data. We could alternatively read
// in the files from disk in the interface layer, and then convert to a map
// in the cpp adapter but that is a bigger faff than just wrapping this
// in a preprocessor flag.
// We'll need something similar down the line when we add Delphi and Python
// interfaces which will have a map over their output data. So we might want
// to revisit this later

template<typename real_type>
using TM1 = Eigen::Tensor<real_type, 1>;

template<typename real_type>
using TM2 = Eigen::Tensor<real_type, 2>;

template<typename real_type>
using TM3 = Eigen::Tensor<real_type, 3>;

template<typename real_type>
using TM4 = Eigen::Tensor<real_type, 4>;

#else

template<typename real_type>
using TM1 = Eigen::TensorMap<Eigen::Tensor<real_type, 1>>;

template<typename real_type>
using TM2 = Eigen::TensorMap<Eigen::Tensor<real_type, 2>>;

template<typename real_type>
using TM3 = Eigen::TensorMap<Eigen::Tensor<real_type, 3>>;

template<typename real_type>
using TM4 = Eigen::TensorMap<Eigen::Tensor<real_type, 4>>;

#endif

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

template<typename real_type>
using T6 = Eigen::Tensor<real_type, 6>;

template<typename real_type>
using T7 = Eigen::Tensor<real_type, 7>;

template<typename real_type, typename std::ptrdiff_t... Dims>
using TFS = Eigen::TensorFixedSize<real_type, Eigen::Sizes<Dims...>>;

}
}

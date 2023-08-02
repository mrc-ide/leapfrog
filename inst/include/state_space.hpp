#pragma once

#include "types.hpp"
#include "model_variants.hpp"

#include <array>

namespace leapfrog {

namespace internal {

template<typename T, std::size_t ... Is>
constexpr std::array<T, sizeof...(Is)>
create_array(T value, std::index_sequence<Is...>) {
  // cast Is to void to remove the warning: unused value
  return {{(static_cast<void>(Is), value)...}};
}

}

template<typename T, std::size_t N>
constexpr std::array<T, N>
create_array(const T &value) {
  return internal::create_array(value, std::make_index_sequence<N>());
}


template<HivAgeStratification S>
struct BaseModelStateSpace;

template<>
struct BaseModelStateSpace<coarse> {
  static constexpr int NS = 2;
  static constexpr int pAG = 81;
  static constexpr int hAG = 9;
  static constexpr int hDS = 7;
  static constexpr int hTS = 3;
  static constexpr std::array<int, 9> hAG_span{2, 3, 5, 5, 5, 5, 5, 5, 31};
};

template<>
struct BaseModelStateSpace<full> {
  static constexpr int NS = 2;
  static constexpr int pAG = 81;
  static constexpr int hAG = 66;
  static constexpr int hDS = 7;
  static constexpr int hTS = 3;
  static constexpr std::array<int, 66> hAG_span = create_array<int, 66>(1);
};

template<bool enabled>
struct ChildModelStateSpace;

template<>
struct ChildModelStateSpace<true> {
  // Number of disease stages within the 1st child age category (0 - 4)
  static constexpr int hc1DS = 6;
  // Number of disease stages within the 2nd child age category (5 - 14)
  static constexpr int hc2DS = 7;
  // Number of transmission types
  static constexpr int hTM = 4;
  // Number of PMTCT types
  static constexpr int hPS = 7;
  // Number of breast feeting age categories
  static constexpr int hBF = 18;
};

template<>
struct ChildModelStateSpace<false> {
};

template<typename ModelVariant>
struct StateSpace {
  static constexpr auto base = BaseModelStateSpace<ModelVariant::stratification>();
  static constexpr auto children = ChildModelStateSpace<ModelVariant::run_child_model>();
};

}

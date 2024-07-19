#pragma once

#include "intermediate_data.hpp"
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
  static constexpr int hc1DS = 7;
  // Number of disease stages within the 2nd child age category (5 - 14)
  static constexpr int hc2DS = 6;
  // last age in the hc1 category
  static constexpr int hc1_ageend = 4;
  // first age in the hc2 category
  static constexpr int hc2_agestart = 5;
  //  Number of age groups in the hc1 category (0-4)
  static constexpr int hc1AG = 5;
  //  Number of age groups in the hc2 category (5-14)
  static constexpr int hc2AG = 10;
  // Age at which the paediatric age group ends (hc1AG + hc2AG)
  static constexpr int hcAG_end = 15;
  // Number of transmission types
  static constexpr int hcTT = 4;
  // Number of PMTCT types
  static constexpr int hPS = 7;
  // Number of PMTCT types with dropout rates
  static constexpr int hPS_dropout = 6;
  // Number of vertical transmission timings (perinatal and breastfeeding)
  static constexpr int hVT = 2;
  // Number of breastfeeding age categories
  static constexpr int hBF = 18;
  // Number of coarse breastfeeding categories
  static constexpr int hBF_coarse = 4;
  // Number of age group possibilities for on ART treatment (0-4,5-9,10-14,0-14)
  static constexpr int hcAG_coarse = 4;
  // Abortion value and indicator of whether it is a percent or number
  static constexpr int hAB_ind = 2;
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

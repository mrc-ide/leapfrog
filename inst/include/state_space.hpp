#pragma once

#include "types.hpp"

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

enum HivAgeStratification {
  full,
  coarse
};

template<HivAgeStratification S>
struct StateSpace;

template<>
struct StateSpace<coarse> {
  static constexpr int num_genders = 2;
  static constexpr int age_groups_pop = 81;
  static constexpr int age_groups_hiv = 9;
  static constexpr int disease_stages = 7;
  static constexpr int treatment_stages = 3;
  static constexpr std::array<int, 9> hiv_age_groups_span{2, 3, 5, 5, 5, 5, 5, 5, 31};
};

template<>
struct StateSpace<full> {
  static constexpr int num_genders = 2;
  static constexpr int age_groups_pop = 81;
  static constexpr int age_groups_hiv = 66;
  static constexpr int disease_stages = 7;
  static constexpr int treatment_stages = 3;
  static constexpr std::array<int, 66> hiv_age_groups_span = create_array<int, 66>(1);
};

template<typename real_type, HivAgeStratification S>
struct Options {
  int hiv_steps;
  double dt;
  const int fertility_first_age_group;
  const int age_groups_fert;
  const int hiv_adult_first_age_group;
  const int adult_incidence_first_age_group;
  const int pAG_INCIDPOP;
  const int time_art_start;
  const int age_groups_hiv_15plus;
  const int scale_cd4_mortality;
  const int hIDX_15PLUS;
  const real_type art_alloc_mxweight;

  Options(double dt,
          int hiv_adult_first_age_group,
          int time_art_start,
          int scale_cd4_mortality,
          real_type art_alloc_mxweight) :
      hiv_steps(hiv_steps),
      dt(dt),
      fertility_first_age_group(15),
      age_groups_fert(35),
      hiv_adult_first_age_group(15),
      adult_incidence_first_age_group(hiv_adult_first_age_group),
      pAG_INCIDPOP(35),
      time_art_start(time_art_start),
      age_groups_hiv_15plus(age_groups_hiv<S>),
      scale_cd4_mortality(scale_cd4_mortality),
      hIDX_15PLUS(0),
      art_alloc_mxweight(art_alloc_mxweight) {}
};

}

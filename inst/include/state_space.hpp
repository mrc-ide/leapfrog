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
}

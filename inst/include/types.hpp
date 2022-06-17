#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

const int MALE = 0;
const int FEMALE = 1;

template <typename real_type>
using TensorMap2 = Eigen::TensorMap<Eigen::Tensor<real_type, 2>>;

template <typename real_type>
using TensorMap3 = Eigen::TensorMap<Eigen::Tensor<real_type, 3>>;

template <typename real_type>
using Tensor2 = Eigen::Tensor<real_type, 2>;

template <typename real_type>
struct Parameters {
  int num_genders;
  int age_groups_pop;             // Default 81 for ages 0 to 80+
  int fertility_first_age_group;  // First index eligible for fertility
  int age_groups_fert;            // Number of ages eligible for fertility

  TensorMap2<real_type> base_pop;
  TensorMap3<real_type> survival;
  TensorMap3<real_type> net_migration;
  TensorMap2<real_type> age_sex_fertility_ratio;
  TensorMap2<real_type> births_sex_prop;
};

template <typename real_type>
struct State {
  Tensor2<real_type> total_population;
  Tensor2<real_type> natural_deaths;
  real_type births;

  State(int age_groups_pop, int num_genders)
      : total_population(age_groups_pop, num_genders),
        natural_deaths(age_groups_pop, num_genders) {}
};

template <typename real_type>
struct WorkingData {
  Tensor2<real_type> migration_rate;

  WorkingData(int age_groups_pop, int num_genders)
      : migration_rate(age_groups_pop, num_genders) {}

  void reset() { migration_rate.setZero(); }
};

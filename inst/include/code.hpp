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

template <typename real_type>
using TensorX2T = Eigen::Tensor<real_type, 2>;

template <typename real_type>
struct Parameters {
  int num_genders;
  int age_groups_pop;             // Default 81 for ages 0 to 80+
  int fertility_first_age_group;  // First index eligible for fertility
  int age_groups_fert;            // Number of ages eligible for fertility

  TensorMapX2T<real_type> base_pop;
  TensorMapX2T<real_type> survival;
  TensorMapX2T<real_type> net_migration;
  TensorMapX1T<real_type> asfr;
  TensorMapX1T<real_type> births_sex_prop;
};

template <typename real_type>
struct State {
  TensorX2T<real_type> total_population;
  TensorX2T<real_type> natural_deaths;
  real_type births;

  State(int age_groups_pop, int num_genders)
      : total_population(age_groups_pop, num_genders),
        natural_deaths(age_groups_pop, num_genders) {}
};

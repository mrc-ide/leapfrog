#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

namespace leapfrog {

template<typename real_type>
using TensorMap1 = Eigen::TensorMap <Eigen::Tensor<real_type, 1>>;

template<typename real_type>
using TensorMap2 = Eigen::TensorMap <Eigen::Tensor<real_type, 2>>;

template<typename real_type>
using TensorMap3 = Eigen::TensorMap <Eigen::Tensor<real_type, 3>>;

template<typename real_type>
using Tensor1 = Eigen::Tensor<real_type, 1>;

template<typename real_type>
using Tensor2 = Eigen::Tensor<real_type, 2>;

template<typename real_type>
using Tensor3 = Eigen::Tensor<real_type, 3>;

template<typename real_type>
using Tensor4 = Eigen::Tensor<real_type, 4>;

template<typename real_type>
struct Parameters {
  int num_genders;
  // Default 81 for ages 0 to 80+
  int age_groups_pop;
  // First index of population eligible for fertility
  int fertility_first_age_group;
  // Number of ages eligible for fertility
  int age_groups_fert;
  // Numer of age groups in HIV population
  int age_groups_hiv;
  // Number of HIV disease stages
  int disease_stages;
  // First index of HIV population to model as adult
  int hiv_adult_first_age_group;
  // Number of HIV ART treatment stages
  int treatment_stages;
  // Time step to start ART treatment
  int time_art_start;
  // Index of the youngest age that is reflected in the adult incidence input
  int adult_incidence_first_age_group;

  int pAG_INCIDPOP;
  // Number of time steps per year in the HIV projection
  int hiv_steps_per_year;
  // Difference in time for each hiv time step in HIV projection
  int dt;

  int scale_cd4_mort;

  // Number of years in each HIV age group
  TensorMap1<int> age_groups_hiv_span;

  // Incidence rate at each time step
  TensorMap1<real_type> incidence_rate;

  TensorMap2<real_type> base_pop;
  TensorMap3<real_type> survival;
  TensorMap3<real_type> net_migration;
  TensorMap2<real_type> age_sex_fertility_ratio;
  TensorMap2<real_type> births_sex_prop;
  TensorMap3<real_type> incidence_relative_risk_age;
  TensorMap1<real_type> incidence_relative_risk_sex;
  TensorMap3<real_type> cd4_mortality;
  TensorMap3<real_type> cd4_progression;
  TensorMap1<int> artcd4elig_idx;
};

template<typename real_type>
struct State {
  Tensor2<real_type> total_population;
  Tensor2<real_type> natural_deaths;
  Tensor2<real_type> hiv_population;
  Tensor2<real_type> hiv_natural_deaths;
  Tensor3<real_type> hiv_strat_adult;
  Tensor4<real_type> art_strat_adult;
  real_type births;
  Tensor3<real_type> aids_deaths_no_art;

  State(int age_groups_pop,
        int num_genders,
        int disease_stages,
        int age_groups_hiv,
        int treatment_stages)
      : total_population(age_groups_pop, num_genders),
        natural_deaths(age_groups_pop, num_genders),
        hiv_population(age_groups_pop, num_genders),
        hiv_natural_deaths(age_groups_pop, num_genders),
        hiv_strat_adult(disease_stages, age_groups_hiv, num_genders),
        art_strat_adult(treatment_stages,
                        disease_stages,
                        age_groups_hiv,
                        num_genders),
        aids_deaths_no_art(disease_stages, age_groups_hiv, num_genders) {}
};

namespace internal {

const int MALE = 0;
const int FEMALE = 1;

template<typename real_type>
struct IntermediateData {
  Tensor2<real_type> migration_rate;
  Tensor2<real_type> hiv_net_migration;
  Tensor2<real_type> hiv_population_coarse_ages;
  Tensor2<real_type> hiv_age_up_prob;
  Tensor2<real_type> hiv_negative_pop;
  Tensor2<real_type> infections_ts;
  Tensor1<real_type> incidence_rate_sex;
  Tensor1<real_type> hiv_neg_aggregate;
  Tensor1<real_type> Xhivn_incagerr;
  Tensor2<real_type> hiv_deaths_age_sex;
  Tensor3<real_type> grad;
  real_type cd4mx_scale;
  real_type artpop_hahm;
  real_type deaths;
  int everARTelig_idx;
  int cd4elig_idx;

  IntermediateData(int age_groups_pop, int age_groups_hiv, int num_genders, int disease_stages)
      : migration_rate(age_groups_pop, num_genders),
        hiv_net_migration(age_groups_pop, num_genders),
        hiv_population_coarse_ages(age_groups_hiv, num_genders),
        hiv_age_up_prob(age_groups_hiv, num_genders),
        hiv_negative_pop(age_groups_pop, num_genders),
        infections_ts(age_groups_pop, num_genders),
        hiv_neg_aggregate(num_genders),
        Xhivn_incagerr(num_genders),
        incidence_rate_sex(num_genders),
        hiv_deaths_age_sex(age_groups_hiv, num_genders),
        grad(disease_stages, age_groups_hiv, num_genders),
        cd4mx_scale(1.0),
        artpop_hahm(0.0),
        deaths(0.0),
        everARTelig_idx(0),
        cd4elig_idx(0) {}

  void reset() {
    migration_rate.setZero();
    hiv_net_migration.setZero();
    hiv_population_coarse_ages.setZero();
    hiv_age_up_prob.setZero();
    hiv_negative_pop.setZero();
    infections_ts.setZero();
    hiv_neg_aggregate.setZero();
    Xhivn_incagerr.setZero();
    incidence_rate_sex.setZero();
    hiv_deaths_age_sex.setZero();
    grad.setZero();
    cd4mx_scale = 1.0;
    deaths = 0.0;
    everARTelig_idx = 0;
    cd4elig_idx = 0;
  }
};

}
}

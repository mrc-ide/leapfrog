#pragma once

#include <model.hpp>
#include <string>

struct Parameters {
  int num_genders = 2;
  int age_groups_pop = 81;             // Default 81 for ages 0 to 80+
  int fertility_first_age_group = 15;  // First index eligible for fertility
  int age_groups_fert = 35;            // Number of ages eligible for fertility

  TensorMapX2cT<real_type> base_pop;
  TensorMapX2cT<real_type> survival;
  TensorMapX2cT<real_type> net_migration;
  TensorMapX1cT<real_type> asfr;
  TensorMapX1cT<real_type> births_sex_prop;
};

struct State {
  TensorMapX2T<real_type> total_population;
  real_type births;
  TensorMapX2T<real_type> natural_deaths;
};

//' The output memory is passed as an argument rather than constructed
//' by the function.
//'
//' Template parameters
//' @param real_type data type variables involved in calculations (typically
//'double' '   but required to be a templated parameter for autodiff
// integration). '
//'
template <typename real_type>
void leapfrog_sim(real_type base_pop,
                  real_type survival,
                  real_type net_migration,
                  real_type asfr,
                  real_type births_sex_group,
                  //
                  // settings
                  std::string model_type,
                  int num_genders,
                  int age_groups_pop,
                  int fertility_first_age_group,
                  int age_groups_fert,
                  int sim_years,
                  //
                  // outputs
                  real_type total_population,
                  real_type births,
                  real_type natural_deaths) {
  const TensorMapX2cT<real_type> base_pop(base_pop, age_groups_pop,
                                          num_genders);
  const TensorMapX2cT<real_type> survival(survival, age_groups_pop + 1,
                                          num_genders);
  const TensorMapX2cT<real_type> net_migration(net_migration, age_groups_pop,
                                               num_genders);
  const TensorMapX1cT<real_type> asfr(asfr, age_groups_fert);
  const TensorMapX1cT<real_type> births_sex_prop(births_sex_group, num_genders);
  Parameters params = Parameters(
      num_genders, age_groups_pop, fertility_first_age_group, age_groups_fert,
      base_pop, survival, net_migration, asfr, births_sex_group);

  // outputs
  TensorMapX2T<real_type> total_population(total_population, age_groups_pop,
                                           num_genders);
  TensorMapX2T<real_type> natural_deaths(natural_deaths, age_groups_pop,
                                         num_genders);
  State initial_state = State(total_population, births, natural_deaths);
  Model model = Model<real_type>(base_pop, survival, net_migration, asfr,
                                 births_sex_group, total_population, births,
                                 natural_deaths);
  model.run_model(sim_years)
}
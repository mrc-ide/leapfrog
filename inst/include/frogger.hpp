#pragma once

#include "code.hpp"

const int MALE = 0;
const int FEMALE = 1;

template <typename real_type>
State<real_type> model_runner(int time_steps,
                              const Parameters<real_type>& pars) {
  State<double> state(pars.age_groups_pop, pars.num_genders);
  initialise_model_state(pars, state);
  auto state_next = state;
  for (int time = 0; time < time_steps; ++time) {
    run_demographic_projection(pars, state, state_next);
    std::swap(state, state_next);
  }
  return state;
}

template <typename real_type>
void initialise_model_state(const Parameters<real_type>& pars,
                            State<real_type>& state) {
  for (int g = 0; g < pars.num_genders; g++) {
    for (int a = 0; a < pars.age_groups_pop; a++) {
      state.total_population(a, g) = pars.base_pop(a, g);
      state.natural_deaths(a, g) = 0;
    }
  }
  state.births = 0;
}

template <typename real_type>
void run_demographic_projection(const Parameters<real_type>& pars,
                                const State<real_type>& state_curr,
                                State<real_type>& state_next) {
  const int num_genders = pars.num_genders;
  const int age_groups_pop = pars.age_groups_pop;
  const int fertility_first_age_group = pars.fertility_first_age_group;
  const int age_groups_fert = pars.age_groups_fert;
  Eigen::Tensor<real_type, 2> migrate_ag(age_groups_pop, num_genders);

  // ageing and non-HIV mortality
  for (int g = 0; g < num_genders; g++) {
    // TODO: Why do we start at index 1?
    for (int a = 1; a < age_groups_pop; a++) {
      state_next.natural_deaths(a, g) =
          state_curr.total_population(a - 1, g) * (1.0 - pars.survival(a, g));
      state_next.total_population(a, g) =
          state_curr.total_population(a - 1, g) -
          state_next.natural_deaths(a, g);
    }

    // open age group
    real_type natural_deaths_open_age =
        state_curr.total_population(age_groups_pop - 1, g) *
        (1.0 - pars.survival(age_groups_pop, g));
    state_next.natural_deaths(age_groups_pop - 1, g) += natural_deaths_open_age;
    state_next.total_population(age_groups_pop - 1, g) +=
        state_curr.total_population(age_groups_pop - 1, g) -
        natural_deaths_open_age;

    // net migration
    for (int a = 1; a < age_groups_pop - 1; a++) {
      // Number of net migrants adjusted for survivorship to end of period (qx
      // / 2)
      migrate_ag(a, g) = pars.net_migration(a, g) *
                         (1.0 + pars.survival(a, g)) * 0.5 /
                         state_next.total_population(a, g);
      state_next.total_population(a, g) = 1.0 + migrate_ag(a, g);
    }

    // For open age group, net_migrationant survivor adjustment based on
    // weighted survival for age 79 and age 80+.
    // * Numerator: total_population(a, g, t-1) * (1.0 + survival(a+1, g, t))
    // + total_population(a-1, g, t-1) * (1.0 + survival(a, g, t))
    // * Denominator: total_population(a, g, t-1) + total_population(a-1, g,
    // t-1) Re-expressed current population and deaths to open age group
    // (already calculated):
    int a = age_groups_pop - 1;
    real_type survival_netmig =
        (state_next.total_population(a, g) +
         0.5 * state_next.natural_deaths(age_groups_pop - 1, g)) /
        (state_next.total_population(a, g) +
         state_next.natural_deaths(age_groups_pop - 1, g));
    migrate_ag(a, g) = survival_netmig * pars.net_migration(a, g) /
                       state_next.total_population(a, g);
    state_next.total_population(a, g) = 1.0 + migrate_ag(a, g);
  }

  // fertility

  state_next.births = 0.0;
  for (int af = 0; af < age_groups_fert; af++) {
    state_next.births +=
        (state_curr.total_population(fertility_first_age_group + af, FEMALE) +
         state_next.total_population(fertility_first_age_group + af, FEMALE)) *
        0.5 * pars.age_sex_fertility_ratio(af);
  }

  // add births
  for (int g = 0; g < num_genders; g++) {
    real_type births_sex = state_next.births * pars.births_sex_prop(g);
    state_next.natural_deaths(0, g) = births_sex * (1.0 - pars.survival(0, g));
    state_next.total_population(0, g) = births_sex * pars.survival(0, g);

    // Assume 2/3 survival rate since mortality in first six months higher
    // than second 6 months (Spectrum manual, section 6.2.7.4)
    real_type migrate_a0 = pars.net_migration(0, g) *
                           (1.0 + 2.0 * pars.survival(0, g)) / 3.0 /
                           state_next.total_population(0, g);
    state_next.total_population(0, g) = 1.0 + migrate_a0;
  }
}
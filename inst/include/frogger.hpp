#pragma once

#include "types.hpp"

const int MALE = 0;
const int FEMALE = 1;

template <typename real_type>
State<real_type> run_model(int time_steps, const Parameters<real_type>& pars) {
  State<double> state(pars.age_groups_pop, pars.num_genders);
  initialise_model_state(pars, state);
  auto state_next = state;
  // Each time step is mid-point of the year
  for (int step = 1; step <= time_steps; ++step) {
    run_demographic_projection(step, pars, state, state_next);
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
void run_demographic_projection(int time_step,
                                const Parameters<real_type>& pars,
                                const State<real_type>& state_curr,
                                State<real_type>& state_next) {
  const int num_genders = pars.num_genders;
  const int age_groups_pop = pars.age_groups_pop;
  const int fertility_first_age_group = pars.fertility_first_age_group;
  const int age_groups_fert = pars.age_groups_fert;
  Tensor2<real_type> migration_rate(age_groups_pop, num_genders);

  // ageing and non-HIV mortality
  for (int g = 0; g < num_genders; ++g) {
    // Start at index 1 as we will add infant (age 0) births and deaths later
    for (int a = 1; a < age_groups_pop; ++a) {
      state_next.natural_deaths(a, g) = state_curr.total_population(a - 1, g) *
                                        (1.0 - pars.survival(a, g, time_step));
      state_next.total_population(a, g) =
          state_curr.total_population(a - 1, g) -
          state_next.natural_deaths(a, g);
    }

    // open age group
    real_type natural_deaths_open_age =
        state_curr.total_population(age_groups_pop - 1, g) *
        (1.0 - pars.survival(age_groups_pop, g, time_step));
    state_next.natural_deaths(age_groups_pop - 1, g) += natural_deaths_open_age;
    state_next.total_population(age_groups_pop - 1, g) +=
        state_curr.total_population(age_groups_pop - 1, g) -
        natural_deaths_open_age;

    // net migration
    // Migration for ages 1, 2, ... 79
    for (int a = 1; a < age_groups_pop - 1; ++a) {
      // Get migration rate, as number of net migrants adjusted for survivorship
      // to end of period. Divide by 2 as (on average) half of deaths will
      // happen before they migrate. Then divide by total pop to get rate.
      migration_rate(a, g) = pars.net_migration(a, g, time_step) *
                             (1.0 + pars.survival(a, g, time_step)) * 0.5 /
                             state_next.total_population(a, g);
      state_next.total_population(a, g) *= 1.0 + migration_rate(a, g);
    }

    // For open age group (age 80+), net migrant survivor adjustment based on
    // weighted survival for age 79 and age 80+.
    // * Numerator: total_population(a, g, t-1) * (1.0 + survival(a+1, g, t))
    // + total_population(a-1, g, t-1) * (1.0 + survival(a, g, t))
    // * Denominator: total_population(a, g, t-1) + total_population(a-1, g,
    // t-1) Re-expressed current population and deaths to open age group
    // (already calculated):
    int a = age_groups_pop - 1;
    real_type survival_netmig =
        (state_next.total_population(a, g) +
         0.5 * state_next.natural_deaths(a, g)) /
        (state_next.total_population(a, g) + state_next.natural_deaths(a, g));
    migration_rate(a, g) = survival_netmig *
                           pars.net_migration(a, g, time_step) /
                           state_next.total_population(a, g);
    state_next.total_population(a, g) *= 1.0 + migration_rate(a, g);
  }

  // fertility

  state_next.births = 0.0;
  for (int af = 0; af < age_groups_fert; ++af) {
    state_next.births +=
        (state_curr.total_population(fertility_first_age_group + af, FEMALE) +
         state_next.total_population(fertility_first_age_group + af, FEMALE)) *
        0.5 * pars.age_sex_fertility_ratio(af, time_step);
  }

  // add births & infant migration
  for (int g = 0; g < num_genders; ++g) {
    real_type births_sex =
        state_next.births * pars.births_sex_prop(g, time_step);
    state_next.natural_deaths(0, g) =
        births_sex * (1.0 - pars.survival(0, g, time_step));
    state_next.total_population(0, g) =
        births_sex * pars.survival(0, g, time_step);

    // Assume 2/3 survival rate since mortality in first six months higher
    // than second 6 months (Spectrum manual, section 6.2.7.4)
    real_type migration_rate_a0 = pars.net_migration(0, g, time_step) *
                                  (1.0 + 2.0 * pars.survival(0, g, time_step)) /
                                  3.0 / state_next.total_population(0, g);
    state_next.total_population(0, g) *= 1.0 + migration_rate_a0;
  }
}

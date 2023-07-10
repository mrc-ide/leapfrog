#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {

template<HivAgeStratification S, typename real_type>
void run_ageing_and_mortality(int time_step,
                              const Parameters<real_type> &pars,
                              const State<S, real_type> &state_curr,
                              State<S, real_type> &state_next,
                              IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto demog = pars.demography;
  for (int g = 0; g < ss.num_genders; ++g) {
    // Start at index 1 as we will add infant (age 0) births and deaths later
    for (int a = 1; a < ss.age_groups_pop; ++a) {
      state_next.natural_deaths(a, g) = state_curr.total_population(a - 1, g) *
                                        (1.0 - demog.survival(a, g, time_step));
      state_next.total_population(a, g) =
          state_curr.total_population(a - 1, g) -
          state_next.natural_deaths(a, g);
    }

    // open age group
    real_type natural_deaths_open_age =
        state_curr.total_population(ss.age_groups_pop - 1, g) *
        (1.0 - demog.survival(ss.age_groups_pop, g, time_step));
    state_next.natural_deaths(ss.age_groups_pop - 1, g) +=
        natural_deaths_open_age;
    state_next.total_population(ss.age_groups_pop - 1, g) +=
        state_curr.total_population(ss.age_groups_pop - 1, g) -
        natural_deaths_open_age;
  }
}

template<HivAgeStratification S, typename real_type>
void run_migration(int time_step,
                   const Parameters<real_type> &pars,
                   const State<S, real_type> &state_curr,
                   State<S, real_type> &state_next,
                   IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto demog = pars.demography;
  for (int g = 0; g < ss.num_genders; ++g) {
    // Migration for ages 1, 2, ... 79
    for (int a = 1; a < ss.age_groups_pop - 1; ++a) {
      // Get migration rate, as number of net migrants adjusted for survivorship
      // to end of period. Divide by 2 as (on average) half of deaths will
      // happen before they migrate. Then divide by total pop to get rate.
      intermediate.migration_rate(a, g) = demog.net_migration(a, g, time_step) *
                                          (1.0 + demog.survival(a, g, time_step)) *
                                          0.5 / state_next.total_population(a, g);
      state_next.total_population(a, g) *= 1.0 + intermediate.migration_rate(a, g);
    }

    // For open age group (age 80+), net migrant survivor adjustment based on
    // weighted survival for age 79 and age 80+.
    // * Numerator: total_population(a, g, t-1) * (1.0 + survival(a+1, g, t))
    // + total_population(a-1, g, t-1) * (1.0 + survival(a, g, t))
    // * Denominator: total_population(a, g, t-1) + total_population(a-1, g,
    // t-1) Re-expressed current population and deaths to open age group
    // (already calculated):
    int a = ss.age_groups_pop - 1;
    real_type survival_netmig =
        (state_next.total_population(a, g) +
         0.5 * state_next.natural_deaths(a, g)) /
        (state_next.total_population(a, g) + state_next.natural_deaths(a, g));
    intermediate.migration_rate(a, g) = survival_netmig *
                                        demog.net_migration(a, g, time_step) /
                                        state_next.total_population(a, g);
    state_next.total_population(a, g) *= 1.0 + intermediate.migration_rate(a, g);
  }
}

template<HivAgeStratification S, typename real_type>
void run_fertility_and_infant_migration(int time_step,
                                        const Parameters<real_type> &pars,
                                        const State<S, real_type> &state_curr,
                                        State<S, real_type> &state_next,
                                        IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto demog = pars.demography;
  state_next.births = 0.0;
  for (int af = 0; af < pars.options.age_groups_fert; ++af) {
    state_next.births += (state_curr.total_population(
        pars.options.fertility_first_age_group + af, FEMALE) +
                          state_next.total_population(
                              pars.options.fertility_first_age_group + af, FEMALE)) *
                         0.5 * demog.age_sex_fertility_ratio(af, time_step);
  }

  // add births & infant migration
  for (int g = 0; g < ss.num_genders; ++g) {
    real_type births_sex =
        state_next.births * demog.births_sex_prop(g, time_step);
    state_next.natural_deaths(0, g) =
        births_sex * (1.0 - demog.survival(0, g, time_step));
    state_next.total_population(0, g) =
        births_sex * demog.survival(0, g, time_step);

    // Assume 2/3 survival rate since mortality in first six months higher
    // than second 6 months (Spectrum manual, section 6.2.7.4)
    real_type migration_rate_a0 = demog.net_migration(0, g, time_step) *
                                  (1.0 + 2.0 * demog.survival(0, g, time_step)) /
                                  3.0 / state_next.total_population(0, g);
    state_next.total_population(0, g) *= 1.0 + migration_rate_a0;
  }
}

}


template<HivAgeStratification S, typename real_type>
void run_general_pop_demographic_projection(int time_step,
                                            const Parameters<real_type> &pars,
                                            const State<S, real_type> &state_curr,
                                            State<S, real_type> &state_next,
                                            internal::IntermediateData<S, real_type> &intermediate) {
  internal::run_ageing_and_mortality<S>(time_step, pars, state_curr, state_next, intermediate);
  internal::run_migration<S>(time_step, pars, state_curr, state_next, intermediate);
  internal::run_fertility_and_infant_migration<S>(time_step, pars, state_curr, state_next, intermediate);
}

}

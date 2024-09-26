#pragma once

#include "intermediate_data.hpp"
#include "generated/state_types.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void run_ageing_and_mortality(int t,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_dm = pars.base.demography;
  const auto& c_ba = state_curr.base;
  auto& n_ba = state_next.base;

  for (int g = 0; g < ss_b.NS; ++g) {
    // Start at index 1 as we will add infant (age 0) births and deaths later
    for (int a = 1; a < ss_b.pAG; ++a) {
      n_ba.p_total_pop_natural_deaths(a, g) = c_ba.p_total_pop(a - 1, g) * (1.0 - p_dm.survival_probability(a, g, t));
      n_ba.p_total_pop(a, g) = c_ba.p_total_pop(a - 1, g) - n_ba.p_total_pop_natural_deaths(a, g);
    }

    // open age group
    real_type p_total_pop_natural_deaths_open_age = c_ba.p_total_pop(ss_b.pAG - 1, g) *
                                                    (1.0 - p_dm.survival_probability(ss_b.pAG, g, t));
    n_ba.p_total_pop_natural_deaths(ss_b.pAG - 1, g) += p_total_pop_natural_deaths_open_age;
    n_ba.p_total_pop(ss_b.pAG - 1, g) += c_ba.p_total_pop(ss_b.pAG - 1, g) - p_total_pop_natural_deaths_open_age;
  }
}

template<typename ModelVariant, typename real_type>
void run_migration(int t,
                   const Parameters<ModelVariant, real_type> &pars,
                   const State<ModelVariant, real_type> &state_curr,
                   State<ModelVariant, real_type> &state_next,
                   IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_dm = pars.base.demography;
  auto& i_ba = intermediate.base;
  auto& n_ba = state_next.base;

  for (int g = 0; g < ss_b.NS; ++g) {
    // Migration for ages 1, 2, ... 79
    for (int a = 1; a < ss_b.pAG - 1; ++a) {
      // Get migration rate, as number of net migrants adjusted for survivorship
      // to end of period. Divide by 2 as (on average) half of deaths will
      // happen before they migrate. Then divide by total pop to get rate.
      i_ba.migration_rate(a, g) = p_dm.net_migration(a, g, t) *
                                  (1.0 + p_dm.survival_probability(a, g, t)) *
                                  0.5 / n_ba.p_total_pop(a, g);
      n_ba.p_total_pop(a, g) *= 1.0 + i_ba.migration_rate(a, g);
    }

    // For open age group (age 80+), net migrant survivor adjustment based on
    // weighted survival_probability for age 79 and age 80+.
    // * Numerator: p_total_pop(a, g, t-1) * (1.0 + survival_probability(a+1, g, t))
    // + p_total_pop(a-1, g, t-1) * (1.0 + survival_probability(a, g, t))
    // * Denominator: p_total_pop(a, g, t-1) + p_total_pop(a-1, g,
    // t-1) Re-expressed current population and deaths to open age group
    // (already calculated):
    int a = ss_b.pAG - 1;
    real_type survival_probability_netmig = (n_ba.p_total_pop(a, g) +
                                            0.5 * n_ba.p_total_pop_natural_deaths(a, g)) /
                                            (n_ba.p_total_pop(a, g) + n_ba.p_total_pop_natural_deaths(a, g));
    i_ba.migration_rate(a, g) = survival_probability_netmig *
                                p_dm.net_migration(a, g, t) /
                                n_ba.p_total_pop(a, g);
    n_ba.p_total_pop(a, g) *= 1.0 + i_ba.migration_rate(a, g);
  }
}

template<typename ModelVariant, typename real_type>
void run_fertility_and_infant_migration(int t,
                                        const Parameters<ModelVariant, real_type> &pars,
                                        const State<ModelVariant, real_type> &state_curr,
                                        State<ModelVariant, real_type> &state_next,
                                        IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  const auto& p_dm = pars.base.demography;
  const auto& p_op = pars.base.options;
  const auto& c_ba = state_curr.base;
  auto& n_ba = state_next.base;

  n_ba.births = 0.0;
  for (int af = 0; af < p_op.p_fertility_age_groups; ++af) {
    auto total_female_pop_per_age_group = c_ba.p_total_pop(p_op.p_idx_fertility_first + af, FEMALE) +
                                          n_ba.p_total_pop(p_op.p_idx_fertility_first + af, FEMALE);
    n_ba.births += total_female_pop_per_age_group * 0.5 * p_dm.age_specific_fertility_rate(af, t);
  }

  // add births & infant migration
  for (int g = 0; g < ss.NS; ++g) {
    real_type births_sex = n_ba.births * p_dm.births_sex_prop(g, t);
    n_ba.p_total_pop_natural_deaths(0, g) = births_sex * (1.0 - p_dm.survival_probability(0, g, t));
    n_ba.p_total_pop(0, g) = births_sex * p_dm.survival_probability(0, g, t);

    if (p_op.proj_period_int == internal::PROJPERIOD_MIDYEAR) {
    // Assume 2/3 survival_probability rate since mortality in first six months higher
    // than second 6 months (Spectrum manual, section 6.2.7.4)
      real_type migration_rate_a0 = p_dm.net_migration(0, g, t) *
                                    (1.0 + 2.0 * p_dm.survival_probability(0, g, t)) /
                                    3.0 / n_ba.p_total_pop(0, g);
      n_ba.p_total_pop(0, g) *= 1.0 + migration_rate_a0;
    }
  }
}

}


template<typename ModelVariant, typename real_type>
void run_general_pop_demographic_projection(int t,
                                            const Parameters<ModelVariant, real_type> &pars,
                                            const State<ModelVariant, real_type> &state_curr,
                                            State<ModelVariant, real_type> &state_next,
                                            internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  const auto& p_op = pars.base.options;
  internal::run_ageing_and_mortality<ModelVariant>(t, pars, state_curr, state_next, intermediate);
  if (p_op.proj_period_int == internal::PROJPERIOD_MIDYEAR) {
    internal::run_migration<ModelVariant>(t, pars, state_curr, state_next, intermediate);
  }
  internal::run_fertility_and_infant_migration<ModelVariant>(t, pars, state_curr, state_next, intermediate);
}


template<typename ModelVariant, typename real_type>
void run_end_year_migration(int t,
                            const Parameters<ModelVariant, real_type> &pars,
                            const State<ModelVariant, real_type> &state_curr,
                            State<ModelVariant, real_type> &state_next,
                            internal::IntermediateData<ModelVariant, real_type> &intermediate) {

  constexpr auto ss_b = StateSpace<ModelVariant>().base;
  const auto& p_dm = pars.base.demography;
  auto& i_ba = intermediate.base;
  auto& n_ba = state_next.base;

  for (int g = 0; g < ss_b.NS; ++g) {

    // Migration for ages 0, ..., 80+
    for (int a = 0; a < ss_b.pAG; ++a) {

      // Calculate migration rate as number of net migrants divided by total pop.
      i_ba.migration_rate(a, g) = p_dm.net_migration(a, g, t) / n_ba.p_total_pop(a, g);
      n_ba.p_total_pop(a, g) *= 1.0 + i_ba.migration_rate(a, g);
    }
  }
}

}

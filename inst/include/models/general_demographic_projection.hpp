#pragma once

#include "../generated/config_mixer.hpp"

namespace leapfrog {
namespace internal {

template<typename Config>
concept GeneralDemographicProjectionEnabled = RunDemographicProjection<Config>;

template<typename Config>
struct GeneralDemographicProjection {
  GeneralDemographicProjection(...) {};
};

template<GeneralDemographicProjectionEnabled Config>
struct GeneralDemographicProjection<Config> {
  using real_type = typename Config::real_type;
  using ModelVariant = typename Config::ModelVariant;
  using SS = Config::SS;
  using Pars = Config::Pars;
  using State = Config::State;
  using Intermediate = Config::Intermediate;
  using Options = Config::Options;
  using Args = Config::Args;

  // private members of this struct
  private:
  // state space
  static constexpr int NS = SS::NS;
  static constexpr int pAG = SS::pAG;
  static constexpr int FEMALE = SS::FEMALE;
  static constexpr int PROJPERIOD_MIDYEAR = SS::PROJPERIOD_MIDYEAR;

  // function args
  int t;
  const Pars& pars;
  const State& state_curr;
  State& state_next;
  Intermediate& intermediate;
  const Options& opts;

  // only exposing the constructor and some methods
  public:
  GeneralDemographicProjection(Args& args):
    t(args.t),
    pars(args.pars),
    state_curr(args.state_curr),
    state_next(args.state_next),
    intermediate(args.intermediate),
    opts(args.opts)
  {};

  void run_general_pop_demographic_projection() {
    run_ageing_and_mortality();
    if (opts.proj_period_int == PROJPERIOD_MIDYEAR) {
      run_migration();
    }
    run_fertility_and_infant_migration();
  };

  void run_end_year_migration() {
    const auto& p_dp = pars.dp;
    auto& i_dp = intermediate.dp;
    auto& n_dp = state_next.dp;

    for (int g = 0; g < NS; ++g) {
      // Migration for ages 0, ..., 80+
      for (int a = 0; a < pAG; ++a) {
        // Calculate migration rate as number of net migrants divided by total pop.
        i_dp.migration_rate(a, g) = p_dp.net_migration(a, g, t) / n_dp.p_total_pop(a, g);
        n_dp.p_total_pop(a, g) *= 1.0 + i_dp.migration_rate(a, g);
      }
    }
  };


  // private methods that we don't want people to call
  private:
  void run_ageing_and_mortality() {
    const auto& p_dp = pars.dp;
    const auto& c_dp = state_curr.dp;
    auto& n_dp = state_next.dp;

    for (int g = 0; g < NS; ++g) {
      // Start at index 1 as we will add infant (age 0) births and deaths later
      for (int a = 1; a < pAG; ++a) {
        n_dp.p_total_pop_natural_deaths(a, g) = c_dp.p_total_pop(a - 1, g) * (1.0 - p_dp.survival_probability(a, g, t));
        n_dp.p_total_pop(a, g) = c_dp.p_total_pop(a - 1, g) - n_dp.p_total_pop_natural_deaths(a, g);
      }

      // open age group
      real_type p_total_pop_natural_deaths_open_age = c_dp.p_total_pop(pAG - 1, g) *
                                                      (1.0 - p_dp.survival_probability(pAG, g, t));
      n_dp.p_total_pop_natural_deaths(pAG - 1, g) += p_total_pop_natural_deaths_open_age;
      n_dp.p_total_pop(pAG - 1, g) += c_dp.p_total_pop(pAG - 1, g) - p_total_pop_natural_deaths_open_age;
    }
  };

  void run_migration() {
    const auto& p_dp = pars.dp;
    auto& n_dp = state_next.dp;
    auto& i_dp = intermediate.dp;

    for (int g = 0; g < NS; ++g) {
      // Migration for ages 1, 2, ... 79
      for (int a = 1; a < pAG - 1; ++a) {
        // Get migration rate, as number of net migrants adjusted for survivorship
        // to end of period. Divide by 2 as (on average) half of deaths will
        // happen before they migrate. Then divide by total pop to get rate.
        i_dp.migration_rate(a, g) = p_dp.net_migration(a, g, t) *
                                    (1.0 + p_dp.survival_probability(a, g, t)) *
                                    0.5 / n_dp.p_total_pop(a, g);
        n_dp.p_total_pop(a, g) *= 1.0 + i_dp.migration_rate(a, g);
      }

      // For open age group (age 80+), net migrant survivor adjustment based on
      // weighted survival_probability for age 79 and age 80+.
      // * Numerator: p_total_pop(a, g, t-1) * (1.0 + survival_probability(a+1, g, t))
      // + p_total_pop(a-1, g, t-1) * (1.0 + survival_probability(a, g, t))
      // * Denominator: p_total_pop(a, g, t-1) + p_total_pop(a-1, g,
      // t-1) Re-expressed current population and deaths to open age group
      // (already calculated):
      int a = pAG - 1;
      real_type survival_probability_netmig = (n_dp.p_total_pop(a, g) +
                                              0.5 * n_dp.p_total_pop_natural_deaths(a, g)) /
                                              (n_dp.p_total_pop(a, g) + n_dp.p_total_pop_natural_deaths(a, g));
      i_dp.migration_rate(a, g) = survival_probability_netmig *
                                  p_dp.net_migration(a, g, t) /
                                  n_dp.p_total_pop(a, g);
      n_dp.p_total_pop(a, g) *= 1.0 + i_dp.migration_rate(a, g);
    }
  };

  void run_fertility_and_infant_migration() {
    const auto& p_dp = pars.dp;
    const auto& c_dp = state_curr.dp;
    auto& n_dp = state_next.dp;

    n_dp.births = 0.0;
    for (int af = 0; af < opts.p_fertility_age_groups; ++af) {
      auto total_female_pop_per_age_group = c_dp.p_total_pop(opts.p_idx_fertility_first + af, FEMALE) +
                                            n_dp.p_total_pop(opts.p_idx_fertility_first + af, FEMALE);
      n_dp.births += total_female_pop_per_age_group * 0.5 * p_dp.age_specific_fertility_rate(af, t);
    }

    // add births & infant migration
    for (int g = 0; g < NS; ++g) {
      real_type births_sex = n_dp.births * p_dp.births_sex_prop(g, t);
      n_dp.p_total_pop_natural_deaths(0, g) = births_sex * (1.0 - p_dp.survival_probability(0, g, t));
      n_dp.p_total_pop(0, g) = births_sex * p_dp.survival_probability(0, g, t);

      if (opts.proj_period_int == PROJPERIOD_MIDYEAR) {
      // Assume 2/3 survival_probability rate since mortality in first six months higher
      // than second 6 months (Spectrum manual, section 6.2.7.4)
        real_type migration_rate_a0 = p_dp.net_migration(0, g, t) *
                                      (1.0 + 2.0 * p_dp.survival_probability(0, g, t)) /
                                      3.0 / n_dp.p_total_pop(0, g);
        n_dp.p_total_pop(0, g) *= 1.0 + migration_rate_a0;
      }
    }
  };
};

}
}

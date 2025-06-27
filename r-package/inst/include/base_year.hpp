#pragma once

#include "generated/model_variants.hpp"
#include "generated/config_mixer.hpp"
#include "options.hpp"

namespace leapfrog {

template<Language L, typename real_type, internal::MV ModelVariant>
void run_base_year_calculations(
  const typename internal::Config<L, real_type, ModelVariant>::Pars& pars,
  typename internal::Config<L, real_type, ModelVariant>::State& initial_state
) {
  using Cfg = internal::Config<L, real_type, ModelVariant>;
  using SS = Cfg::SS;

  // Initialise base year ouputs
  const auto& p_dp = pars.dp;
  auto& is_dp = initial_state.dp;
  const int t0 = 0;

  is_dp.p_total_pop = p_dp.base_pop;

  // Initialise births and deaths in base year. This involves some
  // approximations because the model does not conduct a full
  // cohort component projection for this step.

  // Births in base year
  is_dp.births = 0.0;
  for (int af = 0; af < SS::p_fertility_age_groups; ++af) {

	  // female population exposed to births: average of current age a population
	  // and current age a+1 population survived backwards one year = average of current
	  // and previous year population. (Does not account for any base year net migration)

	  const auto a = SS::p_idx_fertility_first + af;
	  auto female_fertility_population_a = is_dp.p_total_pop(a, SS::FEMALE) + is_dp.p_total_pop(a+1, SS::FEMALE) / p_dp.survival_probability(a+1, SS::FEMALE, t0);
	  is_dp.births += female_fertility_population_a * 0.5 * p_dp.age_specific_fertility_rate(af, t0);
  }

  // Deaths in base year
  // Note: this calculation would **probably** be more accurate if it survived the
  // population 1 year backwards before calculating deaths (discussed with Rob
  // Glaubius 7 June 2025)
  for (int g = 0; g < SS::NS; ++g) {

	  // (a) age 0 deaths
	  is_dp.p_total_pop_background_deaths(0, g) = is_dp.births * p_dp.births_sex_prop(g, t0) *
	    (1.0 - p_dp.survival_probability(0, g, t0));

	  // (b) age 1 to pAG-1 deaths
	  for (int a = 1; a < SS::pAG; ++a) {
	    is_dp.p_total_pop_background_deaths(a, g) = is_dp.p_total_pop(a-1, g) * (1.0 - p_dp.survival_probability(a, g, t0));
	  }

	  // (c) additional deaths from open-ended age group
	  is_dp.p_total_pop_background_deaths(SS::pAG-1, g) += is_dp.p_total_pop(SS::pAG-1, g) * (1.0 - p_dp.survival_probability(SS::pAG-1, g, t0));
  }
}

} // namespace leapfrog

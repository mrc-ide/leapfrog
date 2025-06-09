#pragma once

#include "generated/model_variants.hpp"
#include "generated/config_mixer.hpp"
#include "models/general_demographic_projection.hpp"
#include "models/hiv_demographic_projection.hpp"
#include "models/adult_hiv_model_simulation.hpp"
#include "models/child_model_simulation.hpp"
#include "options.hpp"

#include <format>

namespace leapfrog {

template<Language L, typename real_type, internal::MV ModelVariant>
struct Leapfrog {
  using Cfg = internal::Config<L, real_type, ModelVariant>;
  using SS = Cfg::SS;
  using Pars = Cfg::Pars;
  using State = Cfg::State;
  using Intermediate = Cfg::Intermediate;
  using OutputState = Cfg::OutputState;
  using Args = Cfg::Args;

  static OutputState run_model(
    const Pars& pars,
    const Options<real_type>& opts,
    const std::vector<int> output_years
  ) {
    int simulation_start_year = opts.proj_start_year;

    State initial_state = {};
    initial_state.reset();
    if constexpr (ModelVariant::run_demographic_projection) {

      // Initialise base year ouputs
      const auto& p_dp = pars.dp;
      auto& is_dp = initial_state.dp;
      int t0 = 0;
      
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

    return run_model_from_state(pars, opts, initial_state, simulation_start_year, output_years);
  };

  static OutputState run_model_from_state(
    const Pars& pars,
    const Options<real_type>& opts,
    const State& initial_state,
    const int simulation_start_year,
    const std::vector<int> output_years
  ) {

    const auto min_output_year = std::min_element(std::begin(output_years),
                                                  std::end(output_years));
    if (*min_output_year != opts.proj_start_year && *min_output_year <= simulation_start_year) {
      throw std::invalid_argument(
        std::format("Cannot output year '{}'. Output years must be later than simulation start year '{}'.",
          *min_output_year, simulation_start_year));
    }
    auto state = initial_state;
    auto state_next = state;
    state_next.reset();

    Intermediate intermediate;
    intermediate.reset();

    OutputState output_state(output_years.size());
    save_state(opts.proj_start_year, state, output_state, output_years);

    // Each time step is mid-point of the year
    for (int step = simulation_start_year - opts.proj_start_year + 1; step < opts.proj_steps; ++step) {
      Args args = { step, pars, state, state_next, intermediate, opts };
      project_year(args);
      save_state(opts.proj_start_year + step, state_next,
                 output_state, output_years);
      std::swap(state, state_next);
      state_next.reset();
      intermediate.reset();
    }
    return output_state;
  };

  static State run_model_single_year(
    const Pars& pars,
    const Options<real_type>& opts,
    const State& initial_state,
    const int simulation_start_year
  ) {
    auto state = initial_state;
    auto state_next = state;
    state_next.reset();

    Intermediate intermediate;
    intermediate.reset();

    Args args = { simulation_start_year - opts.proj_start_year + 1, pars, state, state_next, intermediate, opts };
    project_year(args);

    return args.state_next;
  };

  private:
  static void save_state(
    const int step,
    State& state_next,
    OutputState& output_state,
    const std::vector<int>& output_years
  ) {
    for (size_t i = 0; i < output_years.size(); ++i) {
      if (step == output_years[i]) {
        output_state.save_state(i, state_next);
      }
    }
  };

  static void project_year(Args& args) {
    internal::GeneralDemographicProjection<Cfg> general_dp(args);
    internal::HivDemographicProjection<Cfg> hiv_dp(args);
    internal::AdultHivModelSimulation<Cfg> hiv_sim(args);
    internal::ChildModelSimulation<Cfg> child_sim(args);

    if constexpr (ModelVariant::run_demographic_projection) {
      general_dp.run_general_pop_demographic_projection();

      if constexpr (ModelVariant::run_hiv_simulation) {
        hiv_dp.run_hiv_pop_demographic_projection();
        hiv_sim.run_hiv_model_simulation();
      }

      if constexpr (ModelVariant::run_child_model) {
        child_sim.run_child_model_simulation();
      }

      if (args.opts.proj_period_int == SS::PROJPERIOD_CALENDAR) {
        general_dp.run_end_year_migration();

        if constexpr (ModelVariant::run_hiv_simulation) {
          hiv_dp.run_hiv_pop_end_year_migration();
        }
      }
    }
  };
};

}

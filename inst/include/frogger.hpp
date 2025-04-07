#pragma once

#include "generated/model_variants.hpp"
#include "generated/config_mixer.hpp"
#include "models/general_demographic_projection.hpp"
#include "models/hiv_demographic_projection.hpp"
#include "models/adult_hiv_model_simulation.hpp"
#include "models/child_model_simulation.hpp"

namespace leapfrog {

template<typename real_type, internal::MV ModelVariant>
struct Leapfrog {
  using Cfg = Config<real_type, ModelVariant>;
  using SS = Cfg::SS;
  using Pars = Cfg::Pars;
  using State = Cfg::State;
  using Intermediate = Cfg::Intermediate;
  using OutputState = Cfg::OutputState;
  using Options = Cfg::Options;
  using Args = Cfg::Args;

  static const Options get_opts(
    const int hiv_steps,
    const int t_ART_start,
    const bool is_midyear_projection
  ) {
    const int proj_period = is_midyear_projection
      ? internal::BaseSS::PROJPERIOD_MIDYEAR
      : internal::BaseSS::PROJPERIOD_CALENDAR;
    const Options opts = {
      hiv_steps,
      t_ART_start,
      proj_period
    };
    return opts;
  };

  static OutputState run_model(
    const int time_steps,
    const std::vector<int> save_steps,
    const Pars& pars,
    const Options& opts
  ) {
    auto state = State();
    auto state_next = state;

    // TODO Mantra make run time switch with ability to input initial state
    if constexpr (ModelVariant::run_demographic_projection) {
      // set initial state
      state.dp.p_total_pop = pars.dp.base_pop;
    }

    Intermediate intermediate;
    intermediate.reset();

    OutputState output_state(save_steps.size());
    save_state_if_in_save_step(0, state, output_state, save_steps);

    // Each time step is mid-point of the year
    for (int step = 1; step < time_steps; ++step) {
      Args args = { step, pars, state, state_next, intermediate, opts };
      project_year(args);
      save_state_if_in_save_step(step, state_next, output_state, save_steps);
      std::swap(state, state_next);
      state_next.reset();
      intermediate.reset();
    }
    return output_state;
  };

  private:
  static void save_state_if_in_save_step(
    const int step,
    State& state_next,
    OutputState& output_state,
    const std::vector<int>& save_steps
  ) {
    for (size_t i = 0; i < save_steps.size(); ++i) {
      if (step == save_steps[i]) {
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

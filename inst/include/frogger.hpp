#pragma once

#include "general_demographic_projection.hpp"
#include "hiv_demographic_projection.hpp"
#include "model_simulation.hpp"
#include "child_model_simulation.hpp"
#include "state_saver.hpp"
#include "state_space.hpp"
#include "population_adjustment.hpp"


namespace leapfrog {


template<typename ModelVariant, bool pop_adjust, typename real_type>
OutputState<ModelVariant, real_type> run_model(int time_steps,
                                               std::vector<int> save_steps,
                                               const Parameters<ModelVariant, real_type> &pars) {
  auto state = State<ModelVariant, real_type>(pars);
  auto state_next = state;
  const auto demog = pars.base.demography;


  internal::IntermediateData<ModelVariant, real_type> intermediate(pars.base.options.hAG_15plus);

  intermediate.reset();

  StateSaver<ModelVariant, real_type> state_output(time_steps, save_steps);
  // Save initial state
  state_output.save_state(state, 0);

  // Each time step is mid-point of the year
  for (int step = 1; step <= time_steps; ++step) {
    state_next.reset();

    run_general_pop_demographic_projection<ModelVariant>(step, pars, state, state_next,
                                                         intermediate);
    run_hiv_pop_demographic_projection<ModelVariant>(step, pars, state, state_next,
                                                     intermediate);
    run_hiv_model_simulation<ModelVariant>(step, pars, state, state_next, intermediate);
   if constexpr (pop_adjust) {
       run_base_population_adjustment<ModelVariant>(step, pars, state, state_next, intermediate);
    }
    if constexpr (ModelVariant::run_child_model) {
      run_child_model_simulation<ModelVariant>(step, pars, state, state_next, intermediate);
      if constexpr (pop_adjust) {
        run_child_population_adjustment<ModelVariant>(step, pars, state, state_next, intermediate);
     }
    }
    state_output.save_state(state_next, step);
    std::swap(state, state_next);
    intermediate.reset();
  }
  return state_output.get_full_state();
}

}

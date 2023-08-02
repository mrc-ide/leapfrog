#pragma once

#include "general_demographic_projection.hpp"
#include "hiv_demographic_projection.hpp"
#include "model_simulation.hpp"
#include "child_model_simulation.hpp"
#include "state_saver.hpp"
#include "state_space.hpp"

namespace leapfrog {

namespace internal {

template<HivAgeStratification S, typename real_type>
void initialise_model_state(const Parameters<real_type> &pars,
                            State<S, real_type> &state) {
  constexpr auto ss = StateSpace<S>();
  for (int g = 0; g < ss.NS; ++g) {
    for (int a = 0; a < ss.pAG; ++a) {
      state.p_total_pop(a, g) = pars.demography.base_pop(a, g);
    }
  }
  state.p_total_pop_natural_deaths.setZero();
  state.p_hiv_pop.setZero();
  state.p_hiv_pop_natural_deaths.setZero();
  state.h_hiv_adult.setZero();
  state.h_art_adult.setZero();
  state.births = 0;
  state.h_hiv_deaths_no_art.setZero();
  state.p_infections.setZero();
  state.h_hiv_deaths_art.setZero();
  state.h_art_initiation.setZero();
  state.p_hiv_deaths.setZero();
  state.hc_hiv_pop.setZero();
}

}

template<HivAgeStratification S, typename real_type>
typename StateSaver<S, real_type>::OutputState run_model(int time_steps,
                                                         std::vector<int> save_steps,
                                                         const Parameters<real_type> &pars) {
  State<S, real_type> state;

  internal::initialise_model_state<S>(pars, state);
  auto state_next = state;
  internal::IntermediateData<S, real_type> intermediate(pars.options.hAG_15plus);

  intermediate.reset();

  StateSaver<S, real_type> state_output(time_steps, save_steps);
  // Save initial state
  state_output.save_state(state, 0);

  // Each time step is mid-point of the year
  for (int step = 1; step <= time_steps; ++step) {
    state_next.reset();
    run_general_pop_demographic_projection<S>(step, pars, state, state_next, intermediate);
    run_hiv_pop_demographic_projection<S>(step, pars, state, state_next, intermediate);
    run_hiv_model_simulation<S>(step, pars, state, state_next, intermediate);
    if (pars.options.run_child_model) {
      run_child_model_simulation<S>(step, pars, state, state_next, intermediate);
    }
    state_output.save_state(state_next, step);
    std::swap(state, state_next);
    intermediate.reset();
  }
  return state_output.get_full_state();
}

}

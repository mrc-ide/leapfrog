#pragma once

#include "general_demographic_projection.hpp"
#include "hiv_demographic_projection.hpp"
#include "model_simulation.hpp"
#include "state_saver.hpp"
#include "state_space.hpp"

namespace leapfrog {

namespace internal {

template<typename real_type, HivAgeStratification S>
void initialise_model_state(const Parameters<real_type> &pars,
                            State<real_type, S> &state) {
  constexpr auto ss = StateSpace<S>();
  for (int g = 0; g < ss.num_genders; ++g) {
    for (int a = 0; a < ss.age_groups_pop; ++a) {
      state.total_population(a, g) = pars.base_pop(a, g);
    }
  }
  state.natural_deaths.setZero();
  state.hiv_population.setZero();
  state.hiv_natural_deaths.setZero();
  state.hiv_strat_adult.setZero();
  state.art_strat_adult.setZero();
  state.births = 0;
  state.aids_deaths_no_art.setZero();
  state.infections.setZero();
  state.aids_deaths_art.setZero();
  state.art_initiation.setZero();
  state.hiv_deaths.setZero();
}

}

template<typename real_type, HivAgeStratification S>
typename StateSaver<real_type, S>::OutputState run_model(int time_steps,
                                                         std::vector<int> save_steps,
                                                         const Parameters<real_type> &pars) {
  State<real_type, S> state;

  internal::initialise_model_state<real_type, S>(pars, state);
  auto state_next = state;
  internal::IntermediateData<real_type, S> intermediate(pars.age_groups_hiv_15plus);
  intermediate.reset();

  StateSaver<real_type, S> state_output(time_steps, save_steps);
  // Save initial state
  state_output.save_state(state, 0);

  // Each time step is mid-point of the year
  for (int step = 1; step <= time_steps; ++step) {
    state_next.reset();
    run_general_pop_demographic_projection<real_type, S>(step, pars, state, state_next, intermediate);
    run_hiv_pop_demographic_projection<real_type, S>(step, pars, state, state_next, intermediate);
    run_hiv_model_simulation<real_type, S>(step, pars, state, state_next, intermediate);
    state_output.save_state(state_next, step);
    std::swap(state, state_next);
    intermediate.reset();
  }
  return state_output.get_full_state();
}

}

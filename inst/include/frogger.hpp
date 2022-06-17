#pragma once

#include "demographic_projection.hpp"

template <typename real_type>
State<real_type> run_model(int time_steps, const Parameters<real_type>& pars) {
  State<real_type> state(pars.age_groups_pop, pars.num_genders);
  initialise_model_state(pars, state);
  auto state_next = state;
  WorkingData<real_type> working(pars.age_groups_pop, pars.num_genders);
  // Each time step is mid-point of the year
  for (int step = 1; step <= time_steps; ++step) {
    run_demographic_projection(step, pars, state, state_next, working);
    std::swap(state, state_next);
    working.reset();
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

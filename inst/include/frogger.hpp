#pragma once

#include "general_demographic_projection.hpp"
#include "hiv_demographic_projection.hpp"
#include "model_simulation.hpp"
#include "state_saver.hpp"
#include "state_space.hpp"

namespace leapfrog {

namespace internal {

template<typename real_type, HivAgeStratification S>
void initialise_model_state(const StateSpace<S> &ss,
                            const Parameters<real_type> &pars,
                            State<real_type> &state) {
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
typename StateSaver<real_type>::OutputState run_model(int time_steps, const StateSpace<S> &ss,
                                                      std::vector<int> save_steps,
                                                      const Parameters<real_type> &pars) {
  State<real_type> state(ss.age_groups_pop, ss.num_genders,
                         ss.disease_stages, ss.age_groups_hiv,
                         ss.treatment_stages);

  internal::initialise_model_state(ss, pars, state);
  auto state_next = state;
  internal::IntermediateData<real_type> intermediate(ss.age_groups_pop, ss.age_groups_hiv, ss.num_genders,
                                                     ss.disease_stages, ss.treatment_stages,
                                                     pars.age_groups_hiv_15plus);
  intermediate.reset();

  StateSaver<real_type> state_output(time_steps, save_steps, ss.age_groups_pop, ss.num_genders,
                                     ss.disease_stages, ss.age_groups_hiv,
                                     ss.treatment_stages);
  // Save initial state
  state_output.save_state(state, 0);

  // Each time step is mid-point of the year
  for (int step = 1; step <= time_steps; ++step) {
    state_next.reset();
    run_general_pop_demographic_projection(step, ss, pars, state, state_next, intermediate);
    run_hiv_pop_demographic_projection(step, ss, pars, state, state_next, intermediate);
    run_hiv_model_simulation(step, ss, pars, state, state_next, intermediate);
    state_output.save_state(state_next, step);
    std::swap(state, state_next);
    intermediate.reset();
  }
  return state_output.get_full_state();
}

}
